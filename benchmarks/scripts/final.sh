#!/bin/bash
#SBATCH --job-name=FINAL
#SBATCH --output=finallog_%j.txt
#SBATCH --mem=1GB
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia

julia final.jl $SLURM_JOB_ID $1

date

julia clean.jl

