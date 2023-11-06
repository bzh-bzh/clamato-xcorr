import os
import sys
from dotenv import find_dotenv
sys.path.append(os.path.dirname(find_dotenv()))

from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u
import numpy as np
import yaml
import tqdm

import argparse
from pathlib import Path
import glob
import configparser
import random
import multiprocessing

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
parser.add_argument('survey_name', choices=['MOSDEF', 'zDeep', 'VUDS', 'CLAMATO', '3DHST'])
parser.add_argument('--vega-template-folder', type=Path, default=Path('/global/homes/b/bzh/clamato-xcorr/bias-err-pipeline'))
parser.add_argument('--xcorr-folder', type=Path, default=Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'mock' / 'crosscorr')
parser.add_argument('--output-folder', type=Path, default=Path(constants.BIAS_DIR_BASE) / 'xcorr' / 'vega' / 'input')
args = parser.parse_args()


# Define cosmology
cosmo = constants.COSMOLOGY
zmin = 2.0
zmid = 2.3
comdist_mean = cosmo.comoving_distance(zmid)
comdist_zmin = cosmo.comoving_distance(zmin)
dcomdist_dz = cosmo.inv_efunc(zmid) *2998. # in Mpc/h

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

# Read bin edges
PiBin_fil = os.path.join(constants.XCORR_DIR_BASE, 'bins23_pi_0-30hMpc.txt')
SigBin_fil = os.path.join(constants.XCORR_DIR_BASE, 'bins10_sigma_0-30hMpc.txt')

PiBins0 = ascii.read(PiBin_fil)['pi_edges'].data
SigBins0 = ascii.read(SigBin_fil)['sigma_edges'].data

# vega interprets these as center values.
RP_VALS = (PiBins0[1:] + PiBins0[:-1]) / 2
RT_VALS =  (SigBins0[1:] + SigBins0[:-1]) / 2

RT_GRID, RP_GRID = np.meshgrid(RT_VALS, RP_VALS)

Z_EFF = 2.3
N_PROC = multiprocessing.cpu_count()

survey_name = args.survey_name
output_base = args.output_folder
VEGA_TEMPLATE_DIR = args.vega_template_folder

# Case-insensitive glob: https://stackoverflow.com/a/10886685
def insensitive_glob(pattern):
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either, pattern)))

def gen_vega_config_and_data(xcorr_path,
                             use_raw_covar=True):
    
    suffix = Path(xcorr_path).name.replace('xcorrmock', '').replace('.npy', '')
    assert survey_name in suffix
    
    # Generate picca-format FITS file.
    obs_xcorr = np.load(xcorr_path)
    covar_path = insensitive_glob(os.path.join(constants.XCORR_DIR_BASE, 'mock', 'covar', f"covar_{'raw' if use_raw_covar else ''}mock*_"+survey_name+f"_{constants.DATA_VERSION}.npy"))
    assert len(covar_path) == 1
    covar = np.load(covar_path[0])
    assert obs_xcorr.dtype == np.double and covar.dtype == np.double
    
    data_table = Table([covar, obs_xcorr.flatten(), RT_GRID.flatten(order='F'), RP_GRID.flatten(order='F'), np.full(obs_xcorr.size, Z_EFF)],
                      names=['CO', 'DA', 'RT', 'RP', 'Z'])
    data_hdu = fits.BinTableHDU(data=data_table)
    data_hdu.header['RPMIN'] = min(RP_VALS)
    data_hdu.header['RPMAX'] = max(RP_VALS)
    data_hdu.header['NP'] = len(RP_VALS)
    data_hdu.header['RTMIN'] = min(RT_VALS)
    data_hdu.header['RTMAX'] = max(RT_VALS)
    data_hdu.header['NT'] = len(RT_VALS)
    data_hdu.header['BLINDING'] = 'none'
    hdul = fits.HDUList([fits.PrimaryHDU(), data_hdu])
    hdul_path = output_base / f'data{suffix}.fits'
    hdul.writeto(hdul_path, overwrite=True)
    
    # Copy template .ini files to output base folder, and change filenames.
    tracer_cfg = configparser.ConfigParser()
    tracer_cfg.optionxform=str
    tracer_cfg.read(VEGA_TEMPLATE_DIR / 'bias_qsoxlya_template.ini')
    tracer_cfg['data']['filename'] = str(hdul_path)
    tracer_cfg_path = output_base / f'bias_qsoxlya{suffix}.ini'
    with open(tracer_cfg_path, 'w') as f:
        tracer_cfg.write(f)
    
    main_cfg = configparser.ConfigParser()
    main_cfg.optionxform=str
    main_cfg.read(VEGA_TEMPLATE_DIR / 'bias_main_template.ini')
    main_cfg['data sets']['ini files'] = str(tracer_cfg_path)
    
    # Input fixed parameters here.
    main_cfg['parameters']['drp_QSO'] = str(survey_param_dict[survey_name][1])
    main_cfg['parameters']['sigma_velo_disp_gauss_QSO'] = str(survey_param_dict[survey_name][1])
    #main_cfg['parameters']['beta_QSO'] = main_cfg['parameters']['growth_rate'] / main_cfg['parameters']['bias_QSO']
    
    with open(output_base / f'bias_main{suffix}.ini', 'w') as f:
        main_cfg.write(f)
        
xcorr_base = args.xcorr_folder

arguments = sorted(list(glob.glob(str(xcorr_base / '*.npy'))))
(arguments,) = slurm_split(len(arguments), arguments)[0]
#arguments = reversed(arguments)

# with multiprocessing.Pool(N_PROC) as pool:
#     list(tqdm.tqdm(pool.imap(gen_vega_config_and_data, arguments), total=len(arguments)))

for p in tqdm.tqdm(arguments):
    gen_vega_config_and_data(p)