using Test

@testset "ComponentArrays.jl Integration" begin
  include("componentarrays.jl")
end

@testset "Simulation Core" begin
  include("simulation_core.jl")
end

@testset "Open Operators" begin
  include("operators.jl")
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

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

