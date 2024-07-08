#!/bin/bash
#SBATCH --job-name=FINAL
#SBATCH --output=finallog_%j.txt
#SBATCH --mem=500MB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=grauta@ufl.edu
#SBATCH --time=01:00:00

pwd; hostname; date

module load julia

julia $1/final.jl

date
