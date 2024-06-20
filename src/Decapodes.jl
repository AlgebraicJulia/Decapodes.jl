module Decapodes

using ACSets
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using MLStyle

export
gensim, evalsim, compile, compile_env, default_dec_matrix_generate, default_dec_cu_matrix_generate, default_dec_generate,
CPUBackend, CUDABackend, CPUTarget, CUDATarget,
record_dynamics

function record_dynamics()
  error("Please load Makie.jl to use this function")
end;

include("operators.jl")
include("simulation.jl")

# documentation
include("canon/Canon.jl")

end
