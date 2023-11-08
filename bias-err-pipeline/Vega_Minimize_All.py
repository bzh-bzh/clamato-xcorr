import os
import sys
from dotenv import find_dotenv
sys.path.append(os.path.dirname(find_dotenv()))

import numpy as np
import tqdm
from vega import VegaInterface

import argparse
from pathlib import Path
import glob
import configparser
import random
import multiprocessing
import warnings

import constants


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
parser.add_argument('--input-folder', type=Path, default=Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'vega' / 'input')
parser.add_argument('--output-folder', type=Path, default=Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'vega' / 'output')
args = parser.parse_args()

def minimize_galaxy_bias(config_path):    
    suffix = Path(config_path).name.replace('bias_main', '').replace('.ini', '')
    vega = VegaInterface(config_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vega.minimize()
    bestfit_bias = vega.bestfit.params[0].value
    bestfit_bias_err = vega.bestfit.params[0].error
    full_chi2 = vega.chi2()
    reduced_chi2 = vega.chi2() / (vega.data['qsoxlya'].data_size - 1)
    np.save(args.output_folder / f'bestfit_bias{suffix}.npy', np.array([bestfit_bias, bestfit_bias_err, reduced_chi2]))
    
arguments = sorted(list(glob.glob(str(args.input_folder / 'bias_main*.ini'))))
(arguments,) = slurm_split(len(arguments), arguments)[0]

with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    list(tqdm.tqdm(pool.imap(minimize_galaxy_bias, arguments), total=len(arguments)))