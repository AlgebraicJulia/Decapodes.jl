module ConfigGen

using DrWatson
@quickactivate :benchmarks

include(helpersdir("config_helper.jl"))

using Test
using TOML

@testset "Config file validation" begin
  @test is_supported_arch("cpu")
  @test is_supported_arch("cuda")

  @test !is_supported_arch("wrongarch")

  empty_toml = TOML.parse("")
  @test_throws "Configuration is empty" validate_config(empty_toml)

  bad_toml = TOML.parse("[heat.badarch]\ntest_val = [2]\n")
  @test_throws "is not valid" validate_config(bad_toml)

  empty_sim_toml = TOML.parse("[heat.cpu]\n")
  @test_throws "defined but empty" validate_config(empty_sim_toml)

  good_toml = TOML.parse("[physics.cpu]\ntest_val = [2]\n")
  @test validate_config(good_toml) === nothing
end

@testset "Config file loading" begin
  temp_list = Dict()
  temp_params = [Dict("test" => 1)]
  add_meta_data!(temp_list, temp_params)
  @test temp_list["0"]["fields"] == "test"

  temp_list = Dict()
  temp_params = dict_list(Dict("test" => [1, 2]))
  add_task_data!(temp_list, temp_params)
  @test temp_list["1"]["test"] == 1
  @test temp_list["2"]["test"] == 2

  init_params = Dict("full_test" => ["a", "b"])
  temp_list = process_simulation_config(init_params)
  @test temp_list["0"]["fields"] == "full_test"
  @test temp_list["1"]["full_test"] == "a"
  @test temp_list["2"]["full_test"] == "b"
  @test autoconfig_size(temp_list) == 2

  init_params = Dict("full_test" => ["a", "b", "c", "d"])
  temp_list = process_simulation_config(init_params)
  @test autoconfig_size(temp_list) == 4
  @test get_meta_config_info(temp_list)["fields"] == "full_test"
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
