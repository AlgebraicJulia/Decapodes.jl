using Test

@testset "SummationDecapode Construction" begin
  include("summation.jl")
end

@testset "Coordinates" begin
  include("coordinates.jl")
end

@testset "Mesh Loading" begin
  include("meshes.jl")
end

@testset "ComponentArrays.jl Integration" begin
  include("componentarrays.jl")
end

@testset "Simulation" begin
  include("simulation.jl")
end


