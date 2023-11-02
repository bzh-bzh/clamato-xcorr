import os
import sys
from dotenv import find_dotenv
sys.path.append(os.path.dirname(find_dotenv()))

import numpy as np
import tqdm
import time as time
import multiprocessing
import argparse
import glob
import pathlib
import lyafxcorr_kg as xcorr
import constants
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

def slurm_split(n_obj, *args):
    '''Take a variable number of same-length arrays, and split them according to
    environment variables set by SLURM.

    Args:
        n_obj (int): Length of all arrays in *args.

    Returns tuple:
        split_args: Tuple of *args, but split, indexed by the task index.
        offset: int describing array index that our split starts at.
    '''
    if 'SLURM_ARRAY_TASK_COUNT' in os.environ and 'SLURM_ARRAY_TASK_ID' in os.environ:
        use_job_array = True
        n_tasks_env_var = 'SLURM_ARRAY_TASK_COUNT'
        task_ind_env_var = 'SLURM_ARRAY_TASK_ID'
    else:
        use_job_array = False
        n_tasks_env_var = 'SLURM_NTASKS'
        task_ind_env_var = 'SLURM_PROCID'
    if n_tasks_env_var in os.environ and int(os.environ[n_tasks_env_var]) > 1:
        n_tasks = int(os.environ[n_tasks_env_var])
        if use_job_array:
            # Job arrays can have any arbitrary lower index. Safeguard against accidentally submitting it with 1-N values.
            task_ind = int(os.environ[task_ind_env_var]) - int(os.environ['SLURM_ARRAY_TASK_MIN'])
        else:
            task_ind = int(os.environ[task_ind_env_var])
        assert task_ind < n_tasks and n_tasks <= n_obj
        assert np.all([len(a) == n_obj for a in args])
        
        split_args = [np.array_split(a, n_tasks)[task_ind] for a in args]
        split_first_arg = np.array_split(args[0], n_tasks)
        offset = np.sum([len(split_first_arg[i]) for i in range(task_ind)], dtype=int)

        return split_args, offset
    else:
        return args, 0

parser = argparse.ArgumentParser()
parser.add_argument('--gal-folder', type=pathlib.Path, default=pathlib.Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'mock' / 'gal')
parser.add_argument('--output-folder', type=pathlib.Path, default=pathlib.Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'mock' / 'crosscorr')
parser.add_argument('--lya-file', type=pathlib.Path, 
                    default=pathlib.Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'mock' / 'skewers' / 'pixel_radecz_mock_001.bin')
args = parser.parse_args()

# Define cosmology
cosmo = constants.COSMOLOGY

xcorr_dir = constants.XCORR_DIR_BASE
mockdir = os.path.join(xcorr_dir, 'mock')

# Read in bin edges
PiBin_fil = os.path.join(xcorr_dir, 'bins23_pi_0-30hMpc.txt')
SigBin_fil = os.path.join(xcorr_dir, 'bins10_sigma_0-30hMpc.txt')

PiBins0 = ascii.read(PiBin_fil)
SigBins0 = ascii.read(SigBin_fil)

PiEdges = PiBins0['pi_edges'].data
SigEdges = SigBins0['sigma_edges'].data

# Convert bin boundaries from Mpc/h to Mpc
PiEdges  = PiEdges/(len(PiEdges)*[cosmo.h])
SigEdges = SigEdges/(len(SigEdges)*[cosmo.h])

N_PROC = multiprocessing.cpu_count()

lyapix = xcorr.lyapix(args.lya_file, cosmo=cosmo)
npix = lyapix.npix

def gen_crosscorr(galfil):
    #print("Read in {0:d} Ly-a forest pixels from mock {1:03d}".format(lyapix.npix, 
    #                                                                  ivol))
    gal = ascii.read(galfil, format='ipac')
    
    galfil_suffix = pathlib.Path(galfil).name.replace('cat_galmock', '').replace('.dat', '')

    # Convert to 3D Sky positions
    Coord_all = SkyCoord(ra=gal['ra'], dec=gal['dec'],
                         distance=cosmo.comoving_distance(gal['zspec']))
    
    # Compute Cross-Correlations for all galaxies
    XCorr_all, NoXCorr_all = xcorr.xcorr_gal_lya(Coord_all, lyapix, 
                                                   SigEdges, PiEdges, cosmo=cosmo, 
                                                   verbose=0)
    
    np.save(os.path.join(args.output_folder, f'xcorrmock{galfil_suffix}.npy'), XCorr_all.value)

arguments = sorted(list(glob.glob(str(args.gal_folder / '*.dat'))))
(arguments,) = slurm_split(len(arguments), arguments)[0]
#arguments = reversed(arguments)

with multiprocessing.Pool(N_PROC) as pool:
    list(tqdm.tqdm(pool.imap(gen_crosscorr, arguments), total=len(arguments)))