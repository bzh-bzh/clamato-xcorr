"""
Estimate covariance matrix of galaxy-Lya forest cross-correlations via 
bootstrap. 

Takes an argument that points to an ASCII config file. 

Run as:
> python BootCovar_script.py input.cfg
"""

import os
import sys
import numpy as np
import copy
import lyafxcorr_kg as xcorr
import time
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
import yaml

t0 = time.clock()

# Define cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.31)

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
outbase    = input_cfg['output_file_dir_base']
outsuffix  = input_cfg['output_file_suffix']

# If RANDSEED_IN is not None, use it to initialize random number
# generator.

np.random.seed(seed=randseed_in)

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

# This is the loop over the desired resamples
for ii in range(0,nsamp):
    if (ii % 5) == 0:
        print("Iteration #", ii)
    # Make a copy of the pixels and resample
    LyaPixTmp = copy.deepcopy(LyaPix)
    LyaPixTmp = LyaPixTmp.resample()
    # Resample galaxy positions
    GalCoordTmp = GalCoords[np.random.choice(ngal,ngal,replace=True)]
    XCorrTmp, _ = xcorr.xcorr_gal_lya(GalCoordTmp, LyaPixTmp,SigEdges, PiEdges,
                                      cosmo=cosmo,verbose=0)
    XCorrSamples[:,:,ii] = XCorrTmp

# First, reshape the bootstrap samples so that it has dimensions (nbin, nsamp)
XBootReshaped = XCorrSamples.reshape(-1, XCorrSamples.shape[-1])
Covar = np.cov(XBootReshaped)

outfil_cov = os.path.join(outbase, 'covar_'+outsuffix+'.npy')
outfil_samp = os.path.join(outbase, 'bootsamp_'+outsuffix+'.npy')

np.save(outfil_cov, Covar)
np.save(outfil_samp, XCorrSamples)

t1 = time.clock()

print('Total runtime (s) = {}'.format(t1-t0))

exit()
