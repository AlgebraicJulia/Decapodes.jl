#!/bin/bash
#SBATCH --job-name=FINAL
#SBATCH --output=finallog_%j.txt
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=grauta@ufl.edu
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia 

julia --proj=. final.jl table_logs/$SLURM_JOB_ID $1

date

julia clean.jl

