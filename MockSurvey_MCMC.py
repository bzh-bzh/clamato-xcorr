import os
import glob
import sys
import multiprocessing
import pickle
import time
import signal
import logging

from astropy.io import ascii
import emcee
import corner
import yaml
import numpy as np
from scipy.stats import multivariate_normal
from matplotlib import pyplot as plt

import xcorrmodel
import constants


# https://stackoverflow.com/questions/842557/how-to-prevent-a-block-of-code-from-being-interrupted-by-keyboardinterrupt-in-py
class DelayedKeyboardInterrupt:
    def __enter__(self):
        self.signal_received = False
        self.old_handler = signal.signal(signal.SIGINT, self.handler)
                
    def handler(self, sig, frame):
        self.signal_received = (sig, frame)
        logging.debug('SIGINT received. Delaying KeyboardInterrupt.')
    
    def __exit__(self, type, value, traceback):
        signal.signal(signal.SIGINT, self.old_handler)
        if self.signal_received:
            self.old_handler(*self.signal_received)
            
# Subclass of HDFBackend that delays exit upon SIGINT (CTRL-C) until the HDF file is finished being written to,
# so prevent some file consistency errors on sudden interrupt.
# Of course, this won't stop an exit between a call to grow() and a call to save_step(), so it's not totally safe.
class SafeHDFBackend(emcee.backends.HDFBackend):
    def grow(self, *args, **kwargs):
        with DelayedKeyboardInterrupt():
            super().grow(*args, **kwargs)
    
    def save_step(self, *args, **kwargs):
        with DelayedKeyboardInterrupt():
            super().save_step(*args, **kwargs)


os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
# Turn off OMP parallelism as it will interfere with emcee walker parallelism.
# Actually, one of these options causes scipy to single-thread I think.
# TODO: Figure this out?
# os.environ['OMP_NUM_THREADS'] = '1'
# os.environ['OMP_PROC_BIND'] = 'true'
# os.environ['OMP_PLACES'] = 'threads'

# Constants
REGEN_MODEL_CACHE = False

OBS_DIR = os.path.join(constants.XCORR_DIR_BASE, 'obs')
OBS_SUFFIX = '_globalf_'

MOCK_COVAR_DIR = os.path.join(constants.XCORR_DIR_BASE, 'mock', 'covar')
COVAR_SUFFIX = 'mock'

MODEL_DIR = os.path.join(constants.MODEL_DIR_BASE, 'model_grid_v1')

# Read in bin edges
PiBin_fil = os.path.join(constants.XCORR_DIR_BASE, 'bins23_pi_0-30hMpc.txt')
SigBin_fil = os.path.join(constants.XCORR_DIR_BASE, 'bins10_sigma_0-30hMpc.txt')
PiBins0 = ascii.read(PiBin_fil)
SigBins0 = ascii.read(SigBin_fil)
PiEdges = PiBins0['pi_edges'].data
SigEdges = SigBins0['sigma_edges'].data


assert len(sys.argv) == 2
# Open config file and parse parameters
with open(sys.argv[1], 'r') as f:
    input_cfg = yaml.safe_load(f)

survey_name = input_cfg['survey']
init_theta = np.array([input_cfg['prior_bias'], input_cfg['prior_sigz'], input_cfg['prior_dz']]) # Units for sigz and dz are Mpc/h, not cMpc!
assert np.issubdtype(init_theta.dtype, np.number)
n_dim = len(init_theta)
n_walkers = input_cfg['n_walkers']
n_step = input_cfg['n_step']
seed = input_cfg['rng_seed']
# If resuming from a backend hdf5 file, the last step's random state is re-used, so we don't need to worry about
# stale RNG for repeated short-length runs of this script.
np.random.seed(seed)
n_proc = input_cfg['n_processes']
if n_proc == -1:
    n_proc = multiprocessing.cpu_count()


# Load in model grid, as global.
st = time.time()
model_cache_path = os.path.join(constants.MCMC_DIR_BASE, 'model_cache.pickle')
if os.path.exists(model_cache_path) and not REGEN_MODEL_CACHE:
    with open(model_cache_path, 'rb') as f:
        model_list = pickle.load(f)
else:
    model_list = [xcorrmodel.XCorrModel(f) for f in glob.glob(os.path.join(MODEL_DIR, 'linear_cross_*.txt'))]
    with open(model_cache_path, 'wb') as f:
        pickle.dump(model_list, f)
model_func = xcorrmodel.ModelFunc(model_list)
print(f'Loaded model grid in {time.time() - st} sec.')

xcorr_obs = np.load(os.path.join(OBS_DIR, f'xcorr_{survey_name}{OBS_SUFFIX}{constants.DATA_VERSION}.npy'))
covar_path = glob.glob(os.path.join(MOCK_COVAR_DIR, f'covar_{COVAR_SUFFIX}*_{survey_name}_{constants.DATA_VERSION}.npy'))
assert len(covar_path) == 1
covar = np.load(covar_path[0])
invcov = np.linalg.pinv(covar)


# Defined after we load xcorr_obs and invcov so we don't have massive pickling overhead when multiprocessing
def log_likelihood(theta):
    bias, sigz, dz = theta
    xcorr_diff = (xcorr_obs - model_func.XCorrInterpBin(bias, sigz, dz, SigEdges, PiEdges)).flatten()
    chisq = xcorr_diff.T @ invcov @ xcorr_diff
    return -chisq # Ignore constant factors, since M-H only needs a distribution \propto the true prob

def check_bounds(theta):
    return (model_func.bias_lim[0] <= theta[0] <= model_func.bias_lim[1] # bias
          and model_func.sigz_lim[0] <= theta[1] <= model_func.sigz_lim[1] # sig-z
          and -4 <= theta[2] <= 4)

# Flat
def log_prior(theta):
    if check_bounds(theta): # delta-z
        # Extrapolating; don't let the walkers go outside this area by returning -inf.
        return 0
    else:
        return -np.inf

def log_prob(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return log_likelihood(theta) + lp

# print(f'Bounds of bias are [{model_func.bias_lim[0]}, {model_func.bias_lim[1]}]')
# print(f'Bounds of sigz are [{model_func.sigz_lim[0]}, {model_func.sigz_lim[1]}]')

bias_step = 0.1
# sigz_step = np.mean(np.diff(np.sort(model_func.sigz_arr)))
sigz_step = 0.06 # Just eyeballing the smallest sigz step; hopefully the walkers spread out (that's what the docs say).
dz_step = 8 / 81 # From model fitting grid; arbitrary choice
walker_pos = np.random.normal(loc=init_theta, scale=(bias_step, sigz_step, dz_step), size=(n_walkers, n_dim))
for i in range(len(walker_pos)):
    while not check_bounds(walker_pos[i]):
        walker_pos[i] = np.random.normal(loc=init_theta, scale=(bias_step, sigz_step, dz_step))

chain_file_path = os.path.join(constants.MCMC_DIR_BASE, f'chain_{survey_name}.hdf5')
backend = SafeHDFBackend(chain_file_path)
if not os.path.exists(chain_file_path) or (not backend.iteration):
    print(f'Backend file {chain_file_path} is empty: initializing a new chain.')
    backend.reset(n_walkers, n_dim)
else:
    print(f'Resuming from backend file {chain_file_path}, with length {backend.iteration}.')
    walker_pos = None

with multiprocessing.Pool(n_proc) as p:
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob, pool=p, backend=backend)
    sampler.run_mcmc(walker_pos, n_step - backend.iteration, progress=False)
    # intermediate_step = 100
    # for i in np.arange(0, n_step, intermediate_step):
    #     if i == 0:
    #         sampler.run_mcmc(walker_pos, intermediate_step, progress=True)
    #     else:
    #         sampler.run_mcmc(None, intermediate_step, progress=True)
    #     print(f'# iters: {i + intermediate_step}: AT {sampler.get_autocorr_time(quiet=True)} AF {sampler.acceptance_fraction}')