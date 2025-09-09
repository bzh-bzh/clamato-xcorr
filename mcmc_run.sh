#!/bin/bash
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --nodes=5
#SBATCH -t 2:00:00
#SBATCH -o mcmc_run.out
#SBATCH -e mcmc_run.err
#SBATCH --open-mode=append
#SBATCH --mail-user=benjamin.zhang@berkeley.edu
#SBATCH --mail-type=END,FAIL

module unload python
module load python
source ~/.bashrc
conda activate vega
cd /global/u1/b/bzh/clamato-xcorr
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/3DHST.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/CLAMATO.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/MOSDEF.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/VUDS.yaml &
srun -N 1 -n 1 -c 128 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/zDeep.yaml &
wait
