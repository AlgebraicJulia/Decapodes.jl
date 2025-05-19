#!/bin/bash
#SBATCH --job-name=PRINT
#SBATCH --output=printlog_%A_%a.txt
#SBATCH --mem=16GB
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia

SIMNAME=$1
ARCH=$2
TAG=$3

julia array.jl $SLURM_ARRAY_TASK_ID $SIMNAME $ARCH $TAG

date
