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

# Add a depot path

echo "Precompiling Julia"
julia --proj=$DIR --threads=auto -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
echo "Finished precompiling Julia"

mkdir -p results
rm results/*

BENCH_TRACKER="bench_tracker.txt"
CONFIG_TRACKER="config_tracker.txt"
if [ ! -f $BENCH_TRACKER ]; then
  echo "Creating new tracker file"
  echo -n 0 > $BENCH_TRACKER
fi

BENCH_VAL=$(cat $BENCH_TRACKER)
CONFIG_VAL=$(cat $CONFIG_TRACKER)
if [ $BENCH_VAL != $CONFIG_VAL ]; then
  echo "Config has been regenerated, throwing out tuned parameters"
  rm params/*
  echo -n $CONFIG_VAL > $BENCH_TRACKER
fi

mkdir -p params

CPU_JOBS=$(cat count_cpu.txt)

JOBID=$(sbatch --array=1-$CPU_JOBS%6 --parsable --chdir=$DIR array.sh $1)

echo $JOBID

sbatch --dependency=afterok:$JOBID --chdir=$DIR final.sh $1

