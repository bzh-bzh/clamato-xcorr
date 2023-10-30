#!/bin/bash
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --nodes=5
#SBATCH -t 2:00:00
#SBATCH -o bootcovar_all.out
#SBATCH --open-mode=append

module unload python
module load python
source ~/.bashrc
conda activate clamato-xcorr
cd /global/u1/b/bzh/clamato-xcorr
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 BootCovar_script.py bootstrap_cfg/3DHST.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 BootCovar_script.py bootstrap_cfg/CLAMATO.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 BootCovar_script.py bootstrap_cfg/MOSDEF.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 BootCovar_script.py bootstrap_cfg/VUDS.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 BootCovar_script.py bootstrap_cfg/zDeep.yaml &
wait