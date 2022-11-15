using Test

@testset "Composition" begin
  include(joinpath(@__DIR__, "composition.jl"))
end

@testset "Coordinates" begin
  include(joinpath(@__DIR__, "coordinates.jl"))
end

@testset "Diag2DWD" begin
  include(joinpath(@__DIR__, "diag2dwd.jl"))
end

@testset "Meshes" begin
  include(joinpath(@__DIR__, "meshes.jl"))
end

@testset "MultiScaleArrays" begin
  include(joinpath(@__DIR__, "multiscalearrays.jl"))
end

@testset "Summation" begin
  include(joinpath(@__DIR__, "summation.jl"))
end

@testset "Simulation" begin
  include(joinpath(@__DIR__, "simulation.jl"))
end

@testset "Visualization" begin
  include(joinpath(@__DIR__, "visualization.jl"))
@testset "Rewrite" begin
  include(joinpath(@__DIR__, "rewrite.jl"))
end