module Decapodes

using ACSets
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using MLStyle
using NaNMath

export
gensim, evalsim, compile, compile_env, default_dec_matrix_generate, default_dec_cu_matrix_generate, default_dec_generate,
CPUBackend, CUDABackend, CPUTarget, CPUVectorTarget, CUDATarget,
gen_int, eval_int, gen_split, eval_split

include("operators.jl")
include("simulation.jl")

# documentation
include("canon/Canon.jl")

end
