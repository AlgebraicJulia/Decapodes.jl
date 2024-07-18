#!/bin/bash
#SBATCH --job-name=PRINT
#SBATCH --output=printlog_%A_%a.txt
#SBATCH --mem=16GB
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia

julia array.jl $SLURM_ARRAY_TASK_ID $1 $2

date
