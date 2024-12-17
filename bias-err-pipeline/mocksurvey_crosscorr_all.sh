#!/bin/bash
#SBATCH --constraint cpu
#SBATCH --qos regular
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -o /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/mocksurvey_crosscorr_all.out
#SBATCH -e /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/mocksurvey_crosscorr_all.err
#SBATCH --array=0-99

source ~/.bashrc
cd ~/clamato-xcorr/bias-err-pipeline
conda activate clamato-xcorr
python3 MockSurvey_CrossCorr_All.py
