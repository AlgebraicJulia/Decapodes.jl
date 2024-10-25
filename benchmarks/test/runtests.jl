using Test

@testset "Paths and Files" begin
  include("pathfiles.jl")
end

@testset "Main Config Collection" begin
  include("main_config.jl")
end

@testset "Config Generation" begin
  include("config_gen.jl")
end

@testset "Data Aggregation" begin
  include("data_aggr.jl")
end
