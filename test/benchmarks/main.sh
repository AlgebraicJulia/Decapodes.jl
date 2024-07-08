#!/bin/bash
# TODO: Add a precompile script
echo $PWD

DIR=~/git/Decapodes.jl/test/benchmarks
echo $DIR
cd $DIR

echo $PWD

mkdir -p prints
mkdir -p results
mkdir -p params

JOBID=$(sbatch --array=1-4 --parsable array.sh $DIR)

echo $JOBID

sbatch --dependency=afterok:$JOBID final.sh $DIR
