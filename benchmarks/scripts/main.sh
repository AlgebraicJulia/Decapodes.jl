#!/bin/bash
echo $PWD

if [ $# != 1 ]
then
  echo "Usage: 'simulation name'"
  exit 1
fi

SIMULATION=$1

# Add a depot path
DIR=~/git/Decapodes.jl/benchmarks/scripts

echo $DIR
cd $DIR

echo $PWD

module load julia

julia --threads=auto main.jl $SIMULATION

julia --threads=auto clean.jl
