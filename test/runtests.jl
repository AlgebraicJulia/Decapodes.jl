using Test

@testset "Decapodes" begin
  include("Decapodes.jl")
end

@testset "PetriNets" begin
  include("PetriNets.jl")
end

@testset "Examples" begin
  include("Examples.jl")
end

@testset verbose = true "Decapodes2" begin
  include(joinpath(@__DIR__, "Decapodes2", "runtests.jl"))
end
