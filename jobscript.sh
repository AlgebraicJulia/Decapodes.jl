pwd; hostname; date

echo "Hello!"

module load julia

echo "Running Tests..."
julia --project -e 'using Pkg; Pkg.status(); Pkg.test()'

echo "Building Documentation..."
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.status(); Pkg.instantiate(); include("docs/make.jl")'
