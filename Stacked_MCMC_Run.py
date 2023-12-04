import os
import glob
import sys
import multiprocessing
import pickle
import time
import signal
import logging
import collections
import copy
from pathlib import Path

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
# Turn off OMP parallelism as it will interfere with emcee walker parallelism.
# Actually, one of these options causes scipy to single-thread I think.
# TODO: Figure this out?
# os.environ['OMP_NUM_THREADS'] = '1'
# os.environ['OMP_PROC_BIND'] = 'true'
# os.environ['OMP_PLACES'] = 'threads'

from vega import VegaInterface, utils
from astropy.io import ascii
import emcee
import yaml
import numpy as np
from scipy.stats import multivariate_normal
import scipy.linalg
from matplotlib import pyplot as plt

import constants
from MCMC_Run import DelayedKeyboardInterrupt, SafeHDFBackend, EarlyTermEnsembleSampler, one_dimensional_log_lik


# Constants
BASE_DIR = Path('/global/homes/b/bzh/clamato-xcorr/data/split-mcmc')

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

if __name__ == '__main__':
    assert len(sys.argv) == 2
    # Open config file and parse parameters
    with open(sys.argv[1], 'r') as f:
        input_cfg = yaml.safe_load(f)

    survey_name = input_cfg['survey']
    
    assert input_cfg['mass_bin_suffixes']
    
    # We assume the parameter vector theta is formed from the values of an OrderedDict
    # with keys [bias_QSO_{suffix[0]}, bias_QSO_{suffix[1]}, ..., sigma_velo_disp_gauss_QSO, drp_QSO].
    def set_params_all_interfaces(theta, interface_list):
        sigma_velo_disp_gauss_QSO, drp_QSO = theta[-2:]
        bias_list = theta[:-2]
        assert len(bias_list) == len(interface_list)
        for bias, interface in zip(bias_list, interface_list):
            interface.params['sigma_velo_disp_gauss_QSO'] = sigma_velo_disp_gauss_QSO
            interface.params['drp_QSO'] = drp_QSO
            interface.params['bias_QSO'] = bias
    
    # Initialize best-guess vector and parameter limits.
    init_theta = []
    param_limits = []
    if 'limit_overrides' not in input_cfg:
        input_cfg['limit_overrides'] = {}

    def insert_theta_limits(param_name):
        init_theta.append(input_cfg['initial'][param_name])
        if param_name in input_cfg['limit_overrides']:
            param_limits.append(input_cfg['limit_overrides'][param_name])
        else:
            param_limits.append(PARAM_LIMITS_ALL[param_name])

    # Do biases first.
    for _ in range(len(input_cfg['mass_bin_suffixes'])):
        insert_theta_limits('bias_QSO')
    # Now do dispersion, then offset.
    insert_theta_limits('sigma_velo_disp_gauss_QSO')
    insert_theta_limits('drp_QSO')

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

    data_dir = BASE_DIR / survey_name
    
    interface_list = []
    for i, suffix in enumerate(input_cfg['mass_bin_suffixes']):
        vega = VegaInterface(data_dir / f'main_{survey_name}_{suffix}.ini')
        assert not vega.priors

        # Fixed parameters in config.
        if 'fixed' in input_cfg and input_cfg['fixed']:
            for k, p in input_cfg['fixed'].items():
                assert k not in input_cfg['initial']
                vega.params[k] = p
        # Priors in config.
        # Only set priors for one of the interfaces, since these parameters are shared.
        if 'priors' in input_cfg.keys() and input_cfg['priors'] and i == 0:
            assert set(input_cfg['priors'].keys()) == set(['sigma_velo_disp_gauss_QSO', 'drp_QSO'])
            vega.priors = input_cfg['priors']

        # See if we're working in one dimension. If so, monkey-patch the Vega log-likelihood function.
        if input_cfg['one_dimension']:
            VegaInterface.log_lik = one_dimensional_log_lik
        interface_list.append(vega)

    # Defined after we load the vega interface so we don't have massive pickling overhead when multiprocessing.
    def log_likelihood(theta):
        # Since the vega object is copied into each process, we should be able to get away with mutating the params dict stored in it.
        # Not so if we decide to do threading, for whatever reason!
        set_params_all_interfaces(theta, interface_list)
        # This correctly includes normalization between priors and likelihoods.
        return np.sum([interface.log_lik() for interface in interface_list])

    def check_bounds(theta):
        for p, bound in zip(theta, param_limits):
            if not bound[0] <= p <= bound[1]:
                return False
        return True

    # Flat.
    def log_prior(theta):
        if check_bounds(theta):
            return 0
        else:
            # Extrapolating; don't let the walkers go outside this area by returning -inf.
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