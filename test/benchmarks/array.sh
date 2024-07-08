#!/bin/bash
#SBATCH --job-name=PRINT
#SBATCH --output=printlog_%A_%a.txt
#SBATCH --mem=10gb
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia

julia $1/array.jl $SLURM_ARRAY_TASK_ID

date
