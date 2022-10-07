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

@testset "New Decapodes" begin
  include("diag2dwd.jl")
  include("summation.jl")
end

@testset "MultiScaleArrays" begin
  include("multiscalearrays.jl")
end

@testset "Decapodes2" begin
  include(joinpath(@__DIR__, "Decapodes2", "runtests.jl"))
end