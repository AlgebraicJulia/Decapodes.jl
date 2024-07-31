module PathFiles

using DrWatson
@quickactivate :benchmarks

using Test

@testset "File names" begin
  @test "heat_cpu.toml" == get_configname("heat", "cpu")

  @test "stats_1_cpu.jld2" == get_statsfile_name("1", "cpu")
  @test "benchmarks_3_cuda.json" == get_benchfile_name("3", "cuda")

  @test occursin("stats_5_cuda.jld2", get_statsfile("5", "heat", "cuda"))
  @test occursin("benchmarks_7_cpu.json", get_benchfile("7", "heat", "cpu"))

  @test "heat_cpu.toml" == get_configname("heat", "cpu")
  @test "test_cuda.toml" == get_configname("test", "cuda")
  
  expected_config_dir = joinpath("heat", get_configname("heat", "cuda"))
  @test occursin(expected_config_dir, get_config("heat", "cuda"))
end

end