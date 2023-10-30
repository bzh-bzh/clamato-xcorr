#!/bin/bash
#SBATCH -C haswell
#SBATCH -q flex
#SBATCH --nodes=4
#SBATCH -t 8:00:00
#SBATCH --time-min 00:15:00
#SBATCH -o mcmc_run.out
#SBATCH --open-mode=append
#SBATCH --mail-user=benjamin.zhang@berkeley.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --comment=18:00:00  #desired timelimit
#SBATCH --signal=B:USR1@60
#SBATCH --requeue
#SBATCH --open-mode=append

kill_parent() {
	pkill -u bzh -o python3 --signal SIGINT
}

# specify the command to run to checkpoint your job if any (leave blank if none)
ckpt_command=kill_parent

# requeueing the job if reamining time >0 (do not change the following 3 lines )
. /usr/common/software/variable-time-job/setup.sh
requeue_job func_trap USR1
#

module unload python
module load python
source ~/.bashrc
conda activate vega
cd /global/u1/b/bzh/clamato-xcorr
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/3DHST.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/CLAMATO.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/MOSDEF.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/VUDS.yaml &
srun -N 1 -n 1 -c 64 --cpu_bind=cores python3 MCMC_Run.py mcmc_cfg/zDeep.yaml &
wait
