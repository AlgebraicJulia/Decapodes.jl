#!/bin/bash
#SBATCH --job-name=FINAL
#SBATCH --output=finallog_%j.txt
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=grauta@ufl.edu
#SBATCH --time=01:00:00

pwd; hostname; date

cd $1

mkdir -p table_logs/benchmarks_$SLURM_JOB_ID

module load julia 

julia final.jl table_logs/benchmarks_$SLURM_JOB_ID

date

julia clean.jl

