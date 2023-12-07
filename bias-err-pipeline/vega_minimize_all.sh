#!/bin/bash
#SBATCH --constraint cpu
#SBATCH --qos regular
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -o /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/vega_minimize_all.out
#SBATCH -e /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/vega_minimize_all.err
#SBATCH --array=0-9

source ~/.bashrc
cd ~/clamato-xcorr
conda activate vega
python3 bias-err-pipeline/Vega_Minimize_All.py