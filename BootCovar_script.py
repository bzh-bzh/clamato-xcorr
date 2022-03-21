"""
Estimate covariance matrix of galaxy-Lya forest cross-correlations via 
bootstrap. 

Takes an argument that points to an ASCII config file. 

Run as:
> python BootCovar_script.py input.cfg
"""

import multiprocessing
import multiprocessing.dummy
import os
import sys
import numpy as np
import copy
import lyafxcorr_kg as xcorr
import timeit
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
import yaml
import constants

t0 = timeit.default_timer()

# Define cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.31)

assert len(sys.argv) == 2
# Open config file and parse parameters
with open(sys.argv[1], 'r') as f:
    input_cfg = yaml.safe_load(f)

pixfil     = input_cfg['pixel_file']
galfil     = input_cfg['galaxy_file']
cat_str    = input_cfg['survey']
PiBin_fil  = input_cfg['los_bin_file']
SigBin_fil = input_cfg['trans_bin_file']
nsamp      = input_cfg['n_bootstrap']
randseed_in= input_cfg['rng_seed']
nproc      = input_cfg['n_processes']
outbase    = input_cfg['output_file_dir_base']
# outsuffix  = input_cfg['output_file_suffix']
outsuffix = constants.DATA_VERSION

if nproc == -1:
    nproc = multiprocessing.cpu_count()

# Initialize the seed sequence, which is used to generate good (very different PRNG streams) sub-seeds.
seed_seq = np.random.SeedSequence(entropy=randseed_in)
seed_list = seed_seq.spawn(nproc)

# Read in forest pixels
LyaPix = xcorr.lyapix(pixfil,cosmo=cosmo)
print("Read in {0:d} Ly-a forest pixels".format(LyaPix.npix))
npix = LyaPix.npix

# Read in galaxy positions
gal_table = ascii.read(galfil)
gal = gal_table[gal_table['source'] == cat_str]
print('Read in {0:d} galaxies'.format(len(gal)))
assert gal['ra'].unit == u.degree
assert gal['dec'].unit == u.degree

# Generate 3D sky positions for galaxies
GalCoords = SkyCoord(ra=gal['ra'],
                 dec=gal['dec'],
                 distance=cosmo.comoving_distance(
                     gal['zspec']))

# Read in bin edges
PiBins0 = ascii.read(PiBin_fil)
SigBins0 = ascii.read(SigBin_fil)

PiEdges = PiBins0['pi_edges'].data
SigEdges = SigBins0['sigma_edges'].data

# Convert bin boundaries from Mpc/h to Mpc
PiEdges  = PiEdges/(len(PiEdges)*[cosmo.h])
SigEdges = SigEdges/(len(SigEdges)*[cosmo.h])
PiBound = (min(PiEdges), max(PiEdges))

# Initiate Bootstrap
ngal = len(GalCoords)

# Initialize output array to store all the bootstrap samples
XCorrSamples = np.empty([len(SigEdges)-1, len(PiEdges)-1, nsamp])
#LyaPix, GalCoords, SigEdges, PiEdges, cosmo, ngal, 
def do_resampling(n_sub_partial, seed):
    XCorrSamples_Partial = np.zeros([len(SigEdges)-1, len(PiEdges)-1, n_sub_partial])
    rng = np.random.default_rng(seed)
    for ii in range(n_sub_partial):
        # Make a copy of the pixels and resample
        LyaPixTmp = copy.deepcopy(LyaPix)
        LyaPixTmp.rng = rng
        LyaPixTmp = LyaPixTmp.resample()
        # Resample galaxy positions
        GalCoordTmp = GalCoords[rng.choice(ngal,ngal,replace=True)]
        XCorrTmp, _ = xcorr.xcorr_gal_lya(GalCoordTmp, LyaPixTmp,SigEdges, PiEdges,
                                          cosmo=cosmo,verbose=0)
        XCorrSamples_Partial[:, :, ii] = XCorrTmp
    return XCorrSamples_Partial

n_sub_list = [nsamp // nproc for _ in range(nproc)]
n_sub_list[-1] += nsamp % nproc
assert np.sum(n_sub_list) == nsamp
with multiprocessing.Pool(nproc) as pool:
    XCorrSamples = pool.starmap(do_resampling, zip(n_sub_list, seed_list))
XCorrSamples = np.concatenate(XCorrSamples, axis=-1)

# First, reshape the bootstrap samples so that it has dimensions (nbin, nsamp)
XBootReshaped = XCorrSamples.reshape(-1, XCorrSamples.shape[-1])
Covar = np.cov(XBootReshaped)

outfil_cov = os.path.join(outbase, f'covar_{cat_str}_n{nsamp}_{outsuffix}.npy')
outfil_samp = os.path.join(outbase,  f'bootsamp_{cat_str}_n{nsamp}_{outsuffix}.npy')

# DELETE THIS: regression test against old covariances and samples.
# assert np.allclose(Covar, np.load('/global/homes/b/bzh/clamato-xcorr/data/xcorr/bootstrap/covar_mosdef_n100_testserial_old.npy'))
# assert np.allclose(XCorrSamples, np.load('/global/homes/b/bzh/clamato-xcorr/data/xcorr/bootstrap/bootsamp_mosdef_n100_testserial_old.npy'))

np.save(outfil_cov, Covar)
np.save(outfil_samp, XCorrSamples)

t1 = timeit.default_timer()

print('Total runtime (s) = {}'.format(t1-t0))