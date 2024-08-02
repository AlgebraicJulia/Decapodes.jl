#!/bin/bash
#SBATCH --job-name=FINAL
#SBATCH --output=finallog_%j.txt
#SBATCH --mem=1GB
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia

SIMNAME=$1
ARCH=$2
TAG=$3

julia final.jl $SLURM_JOB_ID $SIMNAME $ARCH $TAG

date

julia clean.jl

