using Test

include("pretty.jl")

@testset "Construction" begin
  include("diag2dwd.jl")
end

@testset "SummationDecapode Construction" begin
  include("summation.jl")
end

@testset "Composition" begin
  include("composition.jl")
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

@testset "Operators" begin
  include("operators.jl")
end

#@testset "Average Rewriting" begin
#  include("rewrite.jl")
#end

@testset "Simulation" begin
  include("simulation.jl")
end

@testset "Visualization" begin
  include("visualization.jl")
end

@testset "SummationDecapode Deconstruction" begin
  include("colanguage.jl")
end
