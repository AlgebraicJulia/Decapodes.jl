module PathFiles

using DrWatson
@quickactivate :benchmarks

using Test

@testset "File names" begin
  @test "heat_cpu_test.toml" == simconfig_name(SimNameData("heat", "cpu", "test"))

  @test "stats_heat_cpu_test_1.jld2" == statsfile_name(SimNameData("heat", "cpu", "test", "1"))
  @test "benchmarks_heat_cuda_test_3.json" == benchfile_name(SimNameData("heat", "cuda", "test", "3"))

  @test occursin("stats_heat_cuda_test_5.jld2", statsfile_path(SimNameData("heat", "cuda", "test", "5")))
  @test occursin("benchmarks_heat_cpu_test2_7.json", benchfile_path(SimNameData("heat", "cpu", "test2", "7")))

  @test "heat_cpu_default.toml" == simconfig_name(SimNameData("heat", "cpu", "default"))
  @test "test_cuda_default.toml" == simconfig_name(SimNameData("test", "cuda", "default"))

  expected_config_dir = joinpath("heat", simconfig_name(SimNameData("heat", "cuda", "test")))
  @test occursin(expected_config_dir, simconfig_path(SimNameData("heat", "cuda", "test")))
end

end
