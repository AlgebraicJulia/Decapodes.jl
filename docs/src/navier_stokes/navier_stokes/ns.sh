#!/bin/bash
#SBATCH --job-name=NS                       # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luke.morris@ufl.edu     # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --cpus-per-task=16                  # Number of cores you can use for threads.
#SBATCH --mem=32gb                          # Job memory request
#SBATCH --time=08:00:00                     # Time limit hrs:min:sec
#SBATCH --output=ns_%j.log                  # Standard output and error log
#SBATCH --export=ALL                        # Give the compute node your environment variables (e.g. JULIA_PKG_PRECOMPILE_AUTO)
pwd; hostname; date

module load julia

mkdir -p "/blue/fairbanksj/fairbanksj/jldepot"
export JULIA_DEPOT_PATH="/blue/fairbanksj/fairbanksj/jldepot"
echo "JULIA_DEPOT_PATH:"
echo "$JULIA_DEPOT_PATH"

echo "Launching script"
date

julia --threads=auto --proj=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --threads=auto --proj=. ./ns.jl

date
echo "Exiting script"

