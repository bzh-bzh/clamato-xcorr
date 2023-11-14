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
            
class EarlyTermEnsembleSampler(emcee.EnsembleSampler):            
    def run_mcmc(self, initial_state, nsteps, **kwargs):
        """
        Iterate :func:`sample` for ``nsteps`` iterations and return the result
        Args:
            initial_state: The initial state or position vector. Can also be
                ``None`` to resume from where :func:``run_mcmc`` left off the
                last time it executed.
            nsteps: The number of steps to run.
        Other parameters are directly passed to :func:`sample`.
        This method returns the most recent result from :func:`sample`.
        """
        if initial_state is None:
            if self._previous_state is None:
                raise ValueError(
                    "Cannot have `initial_state=None` if run_mcmc has never "
                    "been called."
                )
            initial_state = self._previous_state

        results = None
        for results in self.sample(initial_state, iterations=nsteps, **kwargs):
            # Changes here.
            if not self.iteration % 100:
                continue
            max_tau = np.max(sampler.get_autocorr_time(tol=0))
            if self.iteration >= 1000 and 1.5 * max_tau < self.iteration // 50:
                print(f'(n_iter // 50): {self.iteration // 50} is >150% of maximum autocorrelation time: {max_tau:.3}; terminating early.')
                break

        # Store so that the ``initial_state=None`` case will work
        self._previous_state = results
        return results

# Constants
BASE_DIR = Path('/global/homes/b/bzh/clamato-xcorr/data/mcmc')

# Order that these are initialized in is the order of the data vector theta which emcee uses.
PARAM_LIMITS_ALL = {
    'bias_QSO': (0, 10),
    'beta_QSO': (0, 10),
    'sigma_velo_disp_gauss_QSO': (0, 20),
    'drp_QSO': (-10, 10),
    'bias_hcd': (-1, 1),
    'beta_hcd': (0, 50),
    'L0_hcd': (0, 20)
}

assert len(sys.argv) == 2
# Open config file and parse parameters
with open(sys.argv[1], 'r') as f:
    input_cfg = yaml.safe_load(f)
# Special handling for beta_QSO.
#if input_cfg['initial']['beta_QSO'] is None:
#    input_cfg['initial']['beta_QSO'] = input_cfg['initial']['bias_QSO']**-1

survey_name = input_cfg['survey']
# Initialize best-guess vector and parameter limits.
init_theta = []
param_limits = collections.OrderedDict()
for k, p in input_cfg['initial'].items():
    if k == 'beta_QSO' and p is None:
        print('beta_QSO initial value is none; fixing inverse relationship with bias_QSO.')
        continue
    init_theta.append(p)
    param_limits[k] = PARAM_LIMITS_ALL[k]
init_theta = np.array(init_theta)
assert np.issubdtype(init_theta.dtype, np.number)
assert len(param_limits) == len(init_theta)
# Parameter limit overrides.
if 'limit_overrides' in input_cfg and input_cfg['limit_overrides']:
    for param_name, lim in input_cfg['limit_overrides'].items():
        print(f'Overriding limit for {param_name} to be {lim}.')
        param_limits[param_name] = tuple(lim)
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
if 'base_dir' in input_cfg:
    BASE_DIR = Path(input_cfg['base_dir'])

# TODO: this is a totally ugly hack to deal with running non-survey-named .ini files. 
# Probably want to think of a better way later.
data_dir = BASE_DIR / survey_name.split('_')[0]
vega = VegaInterface(data_dir / f'main_{survey_name}.ini')
assert not vega.priors
print(np.sum(vega.data['qsoxlya'].mask))

# Fixed parameters in config.
if 'fixed' in input_cfg and input_cfg['fixed']:
    for k, p in input_cfg['fixed'].items():
        assert k not in param_limits
        vega.params[k] = p
# Priors in config.
if 'priors' in input_cfg.keys() and input_cfg['priors']:
    vega.priors = input_cfg['priors']

# Defined after we load the vega interface so we don't have massive pickling overhead when multiprocessing.
def log_likelihood(theta):
    # Since the vega object is copied into each process, we should be able to get away with mutating the params dict stored in it.
    # Not so if we decide to do threading, for whatever reason!
    for p, k in zip(theta, param_limits.keys()):
        vega.params[k] = p
    if not 'beta_QSO' in param_limits.keys():
        vega.params['beta_QSO'] = vega.params['growth_rate'] / vega.params['bias_QSO']
    return -vega.chi2()

def check_bounds(theta):
    for p, bound in zip(theta, param_limits.values()):
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

chain_file_path = data_dir / f'chain_{survey_name}{"_" + input_cfg["chain_file_suffix"] if input_cfg["chain_file_suffix"] else ""}.hdf5'
backend = SafeHDFBackend(chain_file_path)
if not os.path.exists(chain_file_path) or (not backend.iteration):
    print(f'Backend file {chain_file_path} is empty: initializing a new chain.')
    backend.reset(n_walkers, n_dim)
else:
    print(f'Resuming from backend file {chain_file_path}, with length {backend.iteration}.')
    walker_pos = None

with multiprocessing.Pool(n_proc) as p:
    sampler = EarlyTermEnsembleSampler(n_walkers, n_dim, log_prob, pool=p, backend=backend)
    sampler.run_mcmc(walker_pos, n_step - backend.iteration, progress=True)