#!/bin/bash
#SBATCH --constraint cpu
#SBATCH --qos shared
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 1G
#SBATCH -t 01:30:00
#SBATCH -o /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/vega_gen_all.out
#SBATCH -e /global/u1/b/bzh/clamato-xcorr/data/bias-err/papermill/vega_gen_all.err
#SBATCH --array=0-39

source ~/.bashrc
cd ~/clamato-xcorr/bias-err-pipeline
conda activate clamato-xcorr
python3 Vega_Gen_All.py CLAMATO --zero-dispersion true
