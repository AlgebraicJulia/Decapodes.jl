using Test

@info "Executing tests with $(Threads.nthreads()) threads."

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
  # Get the short error description instead of full stacktrace
  error_msg = if isdefined(CUDA, :_initialization_error) && CUDA._initialization_error !== nothing
    CUDA._initialization_error
  else
    "unknown reason"
  end
  @info "CUDA tests were not run, since CUDA.functional() is false." reason=error_msg
end

@testset "Code Quality (Aqua.jl)" begin
  include("aqua.jl")
end

