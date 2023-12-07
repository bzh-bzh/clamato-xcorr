#!/bin/bash
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --nodes=4
#SBATCH -t 4:00:00
#SBATCH -o split_mcmc_run.out
#SBATCH -e split_mcmc_run.err
#SBATCH --open-mode=append
#SBATCH --mail-user=benjamin.zhang@berkeley.edu
#SBATCH --mail-type=END,FAIL

module unload python
module load python
source ~/.bashrc
conda activate vega
cd /global/u1/b/bzh/clamato-xcorr
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 Split_MCMC_Run.py mcmc_cfg/CLAMATO_split.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 Split_MCMC_Run.py mcmc_cfg/MOSDEF_split.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 Split_MCMC_Run.py mcmc_cfg/VUDS_split.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 Split_MCMC_Run.py mcmc_cfg/zDeep_split.yaml &
wait