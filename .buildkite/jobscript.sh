#!/bin/bash

pwd; hostname; date

module load julia

echo "Running Tests..."
julia --project -t 16 -e 'using Pkg; Pkg.status(); Pkg.test()'

echo "Building Documentation..."
julia -t 16 -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.status(); Pkg.up(); include("docs/make.jl")'
