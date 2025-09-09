``#!/bin/bash

module unload python
module load python
source ~/.bashrc
conda activate papermill
#papermill -k clamato-xcorr mocksurvey_lyagen.ipynb out_lya.ipynb
papermill -k clamato-xcorr mocksurvey_galgen.ipynb out_gal.ipynb
conda activate clamato-xcorr
python3 MockSurvey_CrossCorr.py
conda activate papermill
papermill -k clamato-xcorr mocksurvey_covar.ipynb /dev/null
papermill -k clamato-xcorr mcmc_vega_gen.ipynb /dev/null
