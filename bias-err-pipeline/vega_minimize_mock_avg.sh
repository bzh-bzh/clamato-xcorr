#!/bin/bash
#SBATCH --constraint cpu
#SBATCH --qos debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -o /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/vega_minimize_mock_avg.out
#SBATCH -e /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/vega_minimize_mock_avg.err

source ~/.bashrc
cd ~/clamato-xcorr
conda activate vega
python3 bias-err-pipeline/Vega_Minimize_All.py --input-folder /pscratch/sd/b/bzh/clamato-xcorr/mock/crosscorr-avg --output-folder /pscratch/sd/b/bzh/clamato-xcorr/mock/crosscorr-avg --save-bestfit true