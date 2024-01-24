#!/bin/bash

pwd; hostname; date

module load julia

echo "Running Tests..."
julia -t 32 --project -e 'using Pkg; Pkg.status(); Pkg.test()'

echo "Building Documentation..."
julia -t 32 --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.status(); Pkg.instantiate(); include("docs/make.jl")'
