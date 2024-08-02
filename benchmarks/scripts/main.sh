#!/bin/bash

if [ $# != 3 ]
then
  echo "Usage: 'sim_name' 'architecture' 'tag'"
  exit 1
fi

SIMNAME=$1
ARCH=$2
TAG=$3

# TODO: Need a better way to get the scripts dir path
DIR=scripts
cd $DIR

module load julia

julia --threads=auto main.jl $SIMNAME $ARCH $TAG