using Test

@testset "Composition" begin
  include(joinpath(@__DIR__, "composition.jl"))
end

@testset "Rewrite" begin
  include(joinpath(@__DIR__, "rewrite.jl"))
end