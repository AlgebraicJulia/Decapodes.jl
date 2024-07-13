#!/bin/bash
echo $PWD

if [ $# != 1 ]
then
  echo "Usage: 'simulation name'"
  exit 1
fi

DIR=~/git/Decapodes.jl/test/benchmarks

echo $DIR
cd $DIR/$1

echo $PWD

module load julia

# Add a depot path

echo "Precompiling Julia"
julia --proj=$DIR --threads=auto -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
echo "Finished precompiling Julia"

echo "Removing previous results"
mkdir -p results
rm results/*

mkdir -p params

TRACKER="dirty.txt"
if [ -f $TRACKER ]; then
  echo "Config has been regenerated, throwing out tuned parameters"
  rm params/*
  rm $TRACKER
fi

CPU_JOBS=$(cat count_cpu.txt)
echo "Running $CPU_JOBS CPU tasks"
JOBID=$(sbatch --array=1-$CPU_JOBS%6 --parsable --chdir=$DIR array.sh $1)

echo $JOBID

sbatch --dependency=afterok:$JOBID --chdir=$DIR final.sh $1

