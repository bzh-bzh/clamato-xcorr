#!/bin/bash
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --nodes=4
#SBATCH -t 1:30:00
#SBATCH -o mcmc_run_oned.out
#SBATCH -e mcmc_run_oned.err
#SBATCH --open-mode=append
#SBATCH --mail-user=benjamin.zhang@berkeley.edu
#SBATCH --mail-type=END,FAIL

module unload python
module load python
source ~/.bashrc
conda activate vega
cd /global/u1/b/bzh/clamato-xcorr
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/3DHST_oned.yaml &
# srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/CLAMATO_oned.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/MOSDEF_oned.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/VUDS_oned.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/zDeep_oned.yaml &
wait