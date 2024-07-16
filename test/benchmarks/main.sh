#!/bin/bash
echo $PWD

if [ $# != 1 ]
then
  echo "Usage: 'simulation name'"
  exit 1
fi

SIMULATION=$1

DIR=~/git/Decapodes.jl/test/benchmarks
CPU_COUNTER=count_cpu.txt
CUDA_COUNTER=count_cuda.txt

echo $DIR
cd $DIR/$SIMULATION

echo $PWD

if [ ! -f $CPU_COUNTER ] || [ ! -f $CUDA_COUNTER ]; then
  echo "Error: either cpu or cuda benchmarks are not properly configured for $SIMULATION"
  exit 1
fi

module load julia

# Add a depot path

echo "Precompiling Julia"
julia --proj=$DIR --threads=auto -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
echo "Finished precompiling Julia"

echo "Removing previous results"
mkdir -p results
rm --force results/*

mkdir -p params

TRACKER="clean.txt"
if [ ! -f $TRACKER ]; then
  echo "Config has been regenerated, throwing out tuned parameters"
  rm --force params/*
  touch $TRACKER
fi

CPU_JOBS=$(cat $CPU_COUNTER)
echo "Running $CPU_JOBS CPU tasks"
CPU_JOBID=$(sbatch --array=1-$CPU_JOBS%4 --parsable --chdir=$DIR --partition=hpg-milan array.sh $SIMULATION "cpu")
echo "CPU Job ID is: $CPU_JOBID"

CUDA_JOBS=$(cat $CUDA_COUNTER)
echo "Running $CUDA_JOBS CUDA tasks"
CUDA_JOBID=$(sbatch --array=1-$CUDA_JOBS%4 --parsable --chdir=$DIR --partition=gpu --gres=gpu:a100:1 array.sh $SIMULATION "cuda")
echo "CUDA Job ID is: $CUDA_JOBID"

sbatch --dependency=afterok:$CPU_JOBID,$CUDA_JOBID --chdir=$DIR final.sh $SIMULATION

