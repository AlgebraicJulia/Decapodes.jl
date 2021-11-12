using Pkg
Pkg.activate(".")
Pkg.add(url="https://github.com/AlgebraicJulia/CombinatorialSpaces.jl")
Pkg.add(url="https://github.com/bosonbaas/Catlab.jl", rev="rem_fix")
Pkg.develop(url="../..")
Pkg.instantiate()
