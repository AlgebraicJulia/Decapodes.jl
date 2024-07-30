module ConfigGen

using DrWatson
@quickactivate :benchmarks

include(helpersdir("config_helper.jl"))

using Test
using TOML

@testset "Config file validation" begin
  @test validate_arch("cpu")
  @test validate_arch("cuda")

  @test !validate_arch("wrongarch")

  empty_toml = TOML.parse("")
  @test_throws "Configuration is empty" validate_config(empty_toml)

  bad_toml = TOML.parse("[heat.badarch]\ntest_val = [2]\n")
  @test_throws "is not valid" validate_config(bad_toml)

  empty_sim_toml = TOML.parse("[heat.cpu]\n")
  @test_throws "defined but empty" validate_config(empty_sim_toml)

  good_toml = TOML.parse("[physics.cpu]\ntest_val = [2]\n")
  @test validate_config(good_toml) === nothing
end

@testset "Parameter parsing" begin
  test_toml = TOML.parse("[physics.cpu]\ntest_val = [2, 3]\ntest_val2 = [10, 11]")
  @test !isempty(test_toml["physics"]["cpu"])

  test_data = test_toml["physics"]["cpu"]
  task_data = process_simulation_config(test_data)
  @test length(task_data) == 4 + 1 # Each entry plus meta

  for key in keys(task_data)
    key == "0" && continue
    @test keys(task_data[key]) == Set(["test_val", "test_val2"])
  end
end

end