#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH -o Job/job.%j.out
 



module load gcc/11.3.0
module load gurobi/10.0.3

cd FNLBBD


srun ./main 70_$SLURM_ARRAY_TASK_ID 1 10000 100 6 6 1 1 1 0 0 1 0 3601.0