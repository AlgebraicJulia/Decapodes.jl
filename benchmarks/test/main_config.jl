module MainConfig

using DrWatson
@quickactivate :benchmarks

using Test
using TOML

import benchmarks: main_config_physics_info, physics_config_arch_info

@testset "Gather Nested Config Info" begin
  main_config_info = TOML.parse("""
  [heat.cpu.test]
  args = "test"
  [heat.cpu.test2]
  [heat.cuda.cutest]

  [not_heat.cuda.default]
  arg2 = "test2"
  """)

  not_heat_entry = main_config_physics_info(main_config_info, "not_heat")
  @test length(keys(not_heat_entry)) == 1

  not_heat_cuda = physics_config_arch_info(not_heat_entry, "cuda")
  @test keys(not_heat_cuda) == Set(["default"])


  heat_entry = main_config_physics_info(main_config_info, "heat")
  @test length(keys(heat_entry)) == 2

  heat_cpu = physics_config_arch_info(heat_entry, "cpu")
  @test keys(heat_cpu) == Set(["test", "test2"])

  heat_cuda = physics_config_arch_info(heat_entry, "cuda")
  @test keys(heat_cuda) == Set(["cutest"])
end

@testset "Collecting Config Entries" begin
  simple_toml = TOML.parse("""
  [heat.cpu.test]
  """)

  simple_entries = collect_simsfor_physics(simple_toml, "heat")
  @test length(simple_entries) == 1


  multiple_physics_toml = TOML.parse("""
  [heat.cpu.test2]

  [cold.cuda.cutest]
  """)

  multiphys_heat_entries = collect_simsfor_physics(multiple_physics_toml, "heat")
  @test length(multiphys_heat_entries) == 1

  multiphys_cold_entries = collect_simsfor_physics(multiple_physics_toml, "cold")
  @test length(multiphys_heat_entries) == 1

  @test multiphys_cold_entries != multiphys_heat_entries


  multiple_arches_toml = TOML.parse("""
  [heat.cpu.test2]
  [heat.cuda.cutest2]
  [heat.cuda.final]
  """)

  multiarches_entries = collect_simsfor_physics(multiple_arches_toml, "heat")
  @test length(multiarches_entries) == 3
  @test allunique(multiarches_entries)
end

@testset "Automatic Collection" begin
  main_config_info = TOML.parse("""
  [heat.cpu.a]
  args = "test"
  [heat.cpu.b]
  [heat.cuda.c]

  [not_heat.cuda.a]
  arg2 = "test2"
  """)

  hcpa = SimNameData("heat", "cpu", "a")
  @test has_config_args(main_config_info, hcpa)
  hcpa_args = access_config_args(main_config_info, hcpa)
  @test hcpa_args["args"] == "test"

  @test is_supported_arch("cpu")
  @test is_supported_arch("cuda")

  @test !is_supported_arch("wrongarch")

  no_entries = SimNameData("heat", "cuda", "c")
  @test has_config_args(main_config_info, no_entries) == false
  @test_throws "not found in the main config" access_config_args(main_config_info, no_entries)

  no_tag = SimNameData("heat", "cpu", "badtag")
  @test has_config_args(main_config_info, no_tag) == false
  @test_throws "not found in the main config" access_config_args(main_config_info, no_tag)

  no_arch = SimNameData("not_heat", "cpu", "c")
  @test has_config_args(main_config_info, no_arch) == false
  @test_throws "not found in the main config" access_config_args(main_config_info, no_arch)

  bad_arch = SimNameData("not_heat", "fake", "c")
  @test has_config_args(main_config_info, bad_arch) == false
  @test_throws "not found in the main config" access_config_args(main_config_info, bad_arch)

  no_physics = SimNameData("brussel", "cuda", "c")
  @test has_config_args(main_config_info, no_physics) == false
  @test_throws "not found in the main config" access_config_args(main_config_info, no_physics)

  empty_toml = TOML.parse("")
  empty_snd = SimNameData("", "", "")
  @test has_config_args(empty_toml, empty_snd) == false
  @test_throws "not found in the main configuration" access_config_args(empty_toml, empty_snd)
end

end
