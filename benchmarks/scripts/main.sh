#!/bin/bash

DIR=scripts
cd $DIR

module load julia

if [ $# == 3 ]
then
    SIMNAME=$1
    ARCH=$2
    TAG=$3
    julia --threads=auto main.jl $SIMNAME $ARCH $TAG
elif [ $# == 0 ]
then
    julia --threads=auto main.jl
else
  echo "Usage: 'sim_name' 'architecture' 'tag'"
  exit 1
fi

