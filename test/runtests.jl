using Test

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

# @testset "Coordinates" begin
#   include("coordinates.jl")
# end

@testset "ComponentArrays.jl Integration" begin
  include("componentarrays.jl")
end

@testset "Simulation" begin
  include("simulation.jl")
end

using CUDA
if CUDA.functional()
  @testset "CUDA" begin
    include("cuda_sims.jl")
  end
else
  @info "CUDA tests were not run."
  @info CUDA.functional(true)
end
