import os
import sys

import constants

import numpy as np
os.chdir('/global/homes/b/bzh/tardis-tf/tardis')
sys.path.insert(0,'/global/homes/b/bzh/tardis-tf/tardis')
from tardis import *

INPUT_RSD_DENS_PATH = '/global/homes/b/bzh/projectdir/bzh/BolshoiP/z2.5/BolshoiP_rsd_cic_n298_0094.npy'
OUTPUT_PATH_BASE = os.path.join(constants.XCORR_DIR_BASE, 'mock', 'skewers')

SIM_Z = 0.28850**-1 - 1
SNR_MIN, SNR_MAX = 1.4, 10
# In covariance generation, skewers+gals get subsampled for each subsurvey, so this is >CLAMATO density.
SIGHTLINE_SEP_MPC = 1.0

if __name__ == '__main__':
    np.random.seed(352598789)
    os.makedirs(OUTPUT_PATH_BASE, exist_ok=True)
    print(f'Cosmo z: {SIM_Z}')
    rsd_dens = np.load(INPUT_RSD_DENS_PATH)
    # Convert into overdensity 1 + \delta.
    rsd_dens /= np.mean(rsd_dens)
    
    # Make TARDIS universe
    dummy_init_field = np.ones((1,))
    uni_arg_dict = {
        'z_f': SIM_Z,
        'bias_galaxy': [1.49,0.109],
        'bias_lya': [0.226 ,1.5], 
        'bs': 250,
        'nc': rsd_dens.shape[0],
        'dk': 0.1,
        'initial_field': dummy_init_field
    }
    uni = universe(**uni_arg_dict)
    
    # Convert RSD overdens in FGPAed tau, and sub into the TARDIS universe.
    rsd_tau = uni_arg_dict['bias_lya'][0] * rsd_dens**uni_arg_dict['bias_lya'][1]
    print(f'Mean tau is {np.mean(rsd_tau)}')
    print(f'Mean flux is {np.mean(np.exp(-rsd_tau))}')
    rsd_tau = np.expand_dims(rsd_tau, 0)
    uni.final_tau_RSD = rsd_tau

    # Formula from https://physics.stackexchange.com/questions/120093/random-particles-on-a-grid-effect-of-increasing-density-on-distance-between-the
    spec_smoothing_pix = 1.0 * (uni_arg_dict['nc'] / uni_arg_dict['bs'])
    n_skewers = int(round((4 * uni_arg_dict['bs']**2) / (np.pi * SIGHTLINE_SEP_MPC**2)))
    print(f'# skewers: {n_skewers}')
    mock_clamato = lya_survey(uni, SNR_MIN, SNR_MAX, sm=spec_smoothing_pix, n_skewers=n_skewers, add_ondiag_cont_err=True)
    mock_clamato.select_skewers()
    mock_clamato.save_as_dachshund_format(os.path.join(OUTPUT_PATH_BASE, f'allskewers_{uni_arg_dict["nc"]}cubed_{n_skewers}nskew.bin'))