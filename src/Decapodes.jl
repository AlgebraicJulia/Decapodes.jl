module Decapodes

using ACSets
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using MLStyle

export
gensim, evalsim, compile, compile_env, default_dec_matrix_generate, default_dec_cu_matrix_generate, default_dec_generate,
# CartesianPoint, SpherePoint, r, theta, phi, TangentBasis, θhat, ϕhat,
CPUTarget, CUDATarget


# include("coordinates.jl")
include("operators.jl")
include("simulation.jl")

# documentation
include("canon/Canon.jl")

end
