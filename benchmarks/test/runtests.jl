using Test

@testset "Paths and Files" begin
  include("pathfiles.jl")
end

@testset "Config Generation" begin
  include("config_gen.jl")
end