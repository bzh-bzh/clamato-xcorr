#!/bin/bash
#SBATCH --constraint cpu
#SBATCH --qos shared
#SBATCH -n 1
#SBATCH --cpus-per-task 32
#SBATCH --mem-per-cpu 500M
#SBATCH -t 01:30:00
#SBATCH -o /global/u1/b/bzh/clamato-xcorr/data/xcorr/mock/skewers/dachshund/mocksurvey_lya_dachshund.out
#SBATCH -e /global/u1/b/bzh/clamato-xcorr/data/xcorr/mock/skewers/dachshund/mocksurvey_lya_dachshund.err

source ~/.bashrc
for i in {1..189}
do
    /global/u1/b/bzh/dachshund-openmp/dachshund.exe /global/u1/b/bzh/clamato-xcorr/data/xcorr/mock/skewers/dachshund/dachshund_mock_$i.cfg
done
