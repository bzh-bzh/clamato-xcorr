import os
import glob
import sys
import multiprocessing
import pickle
import time
import signal
import logging
import collections
from pathlib import Path

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
# Turn off OMP parallelism as it will interfere with emcee walker parallelism.
# Actually, one of these options causes scipy to single-thread I think.
# TODO: Figure this out?
# os.environ['OMP_NUM_THREADS'] = '1'
# os.environ['OMP_PROC_BIND'] = 'true'
# os.environ['OMP_PLACES'] = 'threads'

from vega import VegaInterface
from astropy.io import ascii
import emcee
import yaml
import numpy as np
from scipy.stats import multivariate_normal
from matplotlib import pyplot as plt

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

# Constants
BASE_DIR = Path('/global/homes/b/bzh/clamato-xcorr/data/mcmc')

# Order that these are initialized in is the order of the data vector theta which emcee uses.
PARAM_LIMITS = collections.OrderedDict()
PARAM_LIMITS['bias_QSO']         = (0, 10)
#PARAM_LIMITS['beta_QSO']         = (0, 10)
PARAM_LIMITS['par_sigma_smooth'] = (0, 10)
PARAM_LIMITS['drp_QSO']          = (-4, 4)
PARAM_LIMITS['bias_hcd']         = (-0.2, 0)
PARAM_LIMITS['beta_hcd']         = (0, 10)

assert len(sys.argv) == 2
# Open config file and parse parameters
with open(sys.argv[1], 'r') as f:
    input_cfg = yaml.safe_load(f)
# Special handling for beta_QSO.
#if input_cfg['priors']['beta_QSO'] is None:
#    input_cfg['priors']['beta_QSO'] = input_cfg['priors']['bias_QSO']**-1

survey_name = input_cfg['survey']
init_theta = np.array([input_cfg['priors'][k] for k in PARAM_LIMITS.keys() if k != 'beta_QSO']) # TODO: remove this?
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
    
data_dir = BASE_DIR / survey_name
vega = VegaInterface(data_dir / f'main_{survey_name}.ini')
assert np.all(vega.data['qsoxlya'].mask)
assert not vega.priors

# Defined after we load the vega interface so we don't have massive pickling overhead when multiprocessing.
def log_likelihood(theta):
    # Since the vega object is copied into each process, we should be able to get away with mutating the params dict stored in it.
    # Not so if we decide to do threading, for whatever reason!
    for p, k in zip(theta, PARAM_LIMITS.keys()):
        vega.params[k] = p
    vega.params['beta_QSO'] = vega.params['bias_QSO']**-1
    return -vega.chi2()

def check_bounds(theta):
    for p, bound in zip(theta, PARAM_LIMITS.values()):
        if not bound[0] <= p <= bound[1]:
            return False
    return True

# Flat
def log_prior(theta):
    if check_bounds(theta):
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

walker_pos = np.random.normal(loc=init_theta, scale=1e-4, size=(n_walkers, n_dim))
for i in range(len(walker_pos)):
    while not check_bounds(walker_pos[i]):
        walker_pos[i] = np.random.normal(loc=init_theta, scale=1e-4, size=(n_dim,))

chain_file_path = data_dir / f'chain_{survey_name}.hdf5'
backend = SafeHDFBackend(chain_file_path)
if not os.path.exists(chain_file_path) or (not backend.iteration):
    print(f'Backend file {chain_file_path} is empty: initializing a new chain.')
    backend.reset(n_walkers, n_dim)
else:
    print(f'Resuming from backend file {chain_file_path}, with length {backend.iteration}.')
    walker_pos = None

with multiprocessing.Pool(n_proc) as p:
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob, pool=p, backend=backend)
    sampler.run_mcmc(walker_pos, n_step - backend.iteration, progress=True)