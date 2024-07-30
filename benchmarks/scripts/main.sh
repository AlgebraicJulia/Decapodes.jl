#!/bin/bash
echo $PWD

if [ $# != 1 ]
then
  echo "Usage: 'simulation name'"
  exit 1
fi

SIMULATION=$1

# TODO: Need a better way to get the scripts dir path
DIR=scripts

echo $DIR
cd $DIR

echo $PWD

module load julia

julia --threads=auto main.jl $SIMULATION