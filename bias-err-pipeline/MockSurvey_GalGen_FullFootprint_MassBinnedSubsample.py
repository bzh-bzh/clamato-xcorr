import os
import sys
from dotenv import find_dotenv
sys.path.append(os.path.dirname(find_dotenv()))

import numpy as np
import numba
import pandas as pd
import time as time
import re
import itertools
import multiprocessing
import constants
import argparse
import pathlib
import lyafxcorr_kg as xcorr

from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

parser = argparse.ArgumentParser()
parser.add_argument('survey_name', choices=['MOSDEF', 'zDeep', 'VUDS', 'CLAMATO', '3DHST'])
parser.add_argument('--output-folder', type=pathlib.Path, default=pathlib.Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'mock' / 'gal')
parser.add_argument('--n-realizations-per-bin', type=int, default=1000)
parser.add_argument('--logmass-lower-bound', type=float, default=9)
parser.add_argument('--logmass-upper-bound', type=float, default=13.4)
parser.add_argument('--logmass-bin-size', type=float, default=0.1)
args = parser.parse_args()

np.random.seed(452813501)
outdir = args.output_folder
os.makedirs(outdir, exist_ok=True)

# Multiprocessing parameters.
N_PROC = 64

# Initialize the seed sequence, which is used to generate good (very different PRNG streams) sub-seeds.
seed_seq = np.random.SeedSequence(entropy=317487248082105438040721345508662758105)

# Define cosmology
cosmo = constants.COSMOLOGY
zmin = 2.0
zmid = 2.3
comdist_mean = cosmo.comoving_distance(zmid)
comdist_zmin = cosmo.comoving_distance(zmin)
dcomdist_dz = cosmo.inv_efunc(zmid) *2998. # in Mpc/h
sim_boxlen_mpch = 250

dz_to_dmpch = lambda z: z * dcomdist_dz
dkms_to_dmpch = lambda v: ((v * u.km / u.s) / cosmo.H(zmid)).value * cosmo.h * (1 + zmid)

# Tuple of (# galaxies, deltaz, sigz).
survey_param_dict = {
    'MOSDEF': (195, -1.12, 1.5),
    'zDeep': (759, -3.04, 2),
    'VUDS': (469, -2.65, dkms_to_dmpch(200)),
    'CLAMATO': (205, -1.27, 1.5),
    '3DHST': (322, -1, dz_to_dmpch(0.0034 * (1 + zmid)))
}

def read_fofp_file(path):
    """
    Read FoF properties file from the give path.
    Just returns a tuple of masses and positions for now, but can be modified
    easily.
    Note that positions are in box size units.

    """
    halo_cat = pd.read_csv(path, header=0, skiprows=30, delim_whitespace=True)
    halo_cat.rename(columns={'#ID': 'ID'}, inplace=True)
    print(halo_cat.columns)
    print("{} halos read".format(len(halo_cat)))

    with open(path, 'r') as f:
        for i, l in enumerate(f.readlines()):
            match = re.search(r'#h = ([0-9].[0-9]*)', l)
            if match is not None:
                littleh = float(match.group(1))
                print(f'Found little h to be {littleh}')
                break
            elif i > 200:
                raise RuntimeError
                
    # The mass is initially in Msun (no h). We DO NOT convert to Msun/h by dividing by the h provided in the file.
    # Halo virial mass
    halo_masses = halo_cat.M.to_numpy() #* littleh
    # Stellar mass (truth, including intrinsic spread?)
    stellar_masses_true = halo_cat.SM.to_numpy() #* littleh
    # Stellar mass (including observational + systematic errors)
    stellar_masses_obs = halo_cat.obs_SM.to_numpy() #* littleh
    
    # sigma = np.fromfile(fof_file, dtype="f4", count=num_groups)
    # v_circ = np.fromfile(fof_file, dtype="f4", count=num_groups)
    # min_id = np.fromfile(fof_file, dtype="f4", count=num_groups)
    # v_pot_e = np.fromfile(fof_file, dtype="f4", count=num_groups)
    
    # Apply RSD.
    z_rsd_offset = halo_cat.VZ * ((1 + zmid) / cosmo.H(zmid) * littleh).value
    halo_cat.Z += z_rsd_offset
    
    # In Mpc/h
    positions = np.array([halo_cat.X, halo_cat.Y, halo_cat.Z]).T

    return (halo_masses, stellar_masses_true, stellar_masses_obs, positions)

fofp_path = os.path.join(constants.SIM_DIR_BASE, 'sfr_catalog_0.288498.txt')

masses, _, _, positions = read_fofp_file(fofp_path)

# Generate mock redshift catalogs

# Grab all halos within this footprint
xhalos = positions[:,0]#*256.
yhalos = positions[:,1]#*256.
zhalos = positions[:,2]#*256.

# @numba.njit
def select_halos(logmass_lb, logmass_ub, n_halos_subsamp, delta_z, sigma_z, seed):
    # This should be different for every function call!
    np.random.seed(seed)
    
    mass_mask_indices = np.argwhere(np.logical_and(masses > 10**logmass_lb, masses <= 10**logmass_ub)).flatten()
    
    # Select n_halos_subsamp halos from halos within the mass range.
    halos_here_indices = np.random.choice(mass_mask_indices, 
                                          size=n_halos_subsamp, 
                                          # Sample with replacement since we want bootstrap errors.
                                          replace=True)
    nh_tmp = n_halos_subsamp
    #print('{} halos within footprint'.format(nh_tmp))

    # Positions and stellar masses of all halos within footprint
    xh_tmp = xhalos[halos_here_indices]
    yh_tmp = yhalos[halos_here_indices]
    zh_tmp = zhalos[halos_here_indices]
        
    # Here, we also move half of the halos to the extended, second z-half of the volume
    ind_half_tmp = np.random.choice(np.arange(nh_tmp), 
                                    size=nh_tmp // 2, 
                                    replace=False)
    zh_tmp[ind_half_tmp] = zh_tmp[ind_half_tmp] + sim_boxlen_mpch

    # For z-dimension, also introduce redshift offset and scatter
    zh_tmp += np.random.normal(delta_z, sigma_z, size=len(zh_tmp))

    return (xh_tmp, yh_tmp, zh_tmp)


def gen_gal_realization(survey_name, logmass_lb, logmass_ub, seed_seq_child):
    for i in np.arange(args.n_realizations_per_bin):
        
        # Numba has no support for using a RandomState instance, so we have to reseed every single function call.
        get_new_seed = lambda: seed_seq_child.spawn(1)[0].generate_state(4)[0]
        
        x_all, y_all, z_all = select_halos(logmass_lb, logmass_ub, 
                                                 *survey_param_dict[survey_name], 
                                                 get_new_seed())

        # Convert 3D positions in RA, Dec, and redshift
        ra_all = 180./np.pi * x_all/comdist_mean.value/cosmo.h
        dec_all = 180./np.pi * y_all/comdist_mean.value/cosmo.h
        red_all = zmin + z_all/dcomdist_dz
        # Generate array of survey names.
        n_halos_subsamp = survey_param_dict[survey_name][0]
        source_all = np.asarray([survey_name] * n_halos_subsamp)
        smass_all = np.ones(n_halos_subsamp) * np.nan

        # Output to file
        output_table = Table([ra_all, dec_all, red_all, source_all, smass_all],
                            names=('ra', 'dec', 'zspec', 'source', 'stellar_mass'),
                            dtype=('f8', 'f8', 'f4', 'U', 'f8') )
        output_table['ra'].unit=u.degree
        output_table['dec'].unit = u.degree
        outname = os.path.join(outdir, 'cat_galmock_survey_{}_mass_{}_{}_realization_{}.dat'.format(survey_name, logmass_lb, logmass_ub, i))
        ascii.write(output_table, outname, format='ipac', overwrite=True)

gal_realization_args = []
log_lb_list = np.arange(args.logmass_lower_bound, args.logmass_upper_bound - args.logmass_bin_size, args.logmass_bin_size)
seed_list = seed_seq.spawn(len(log_lb_list))
for log_lb, seed in zip(log_lb_list, seed_list):
    gal_realization_args.append((args.survey_name, log_lb, log_lb + args.logmass_bin_size, seed))

with multiprocessing.Pool(N_PROC) as pool:
    pool.starmap(gen_gal_realization, gal_realization_args)