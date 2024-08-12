module MainConfig

using DrWatson
@quickactivate :benchmarks

using Test
using TOML

include(helpersdir("main_config_helper.jl"))

@testset "Gather Nested Config Info" begin
  main_config_info = TOML.parse("""
  [heat]
  cpu = ["test", "test2"]
  cuda = ["cutest"]

  [not_heat]
  cuda = ["default"]
  """)

  not_heat_entry = main_config_physics_info(main_config_info, "not_heat")
  @test length(keys(not_heat_entry)) == 1

  not_heat_cuda = physics_config_arch_info(not_heat_entry, "cuda")
  @test not_heat_cuda == ["default"]


  heat_entry = main_config_physics_info(main_config_info, "heat")
  @test length(keys(heat_entry)) == 2

  heat_cpu = physics_config_arch_info(heat_entry, "cpu")
  @test heat_cpu == ["test", "test2"]

  heat_cuda = physics_config_arch_info(heat_entry, "cuda")
  @test heat_cuda == ["cutest"]
end

@testset "Collecting Config Entries" begin
  simple_toml = TOML.parse("""
  [heat]
  cpu = ["test"]
  """)

  simple_entries = collect_entriesfor_physics(simple_toml, "heat")
  @test length(simple_entries) == 1


  multiple_physics_toml = TOML.parse("""
  [heat]
  cpu = ["test2"]

  [cold]
  cuda = ["cutest"]
  """)

  multiphys_heat_entries = collect_entriesfor_physics(multiple_physics_toml, "heat")
  @test length(multiphys_heat_entries) == 1

  multiphys_cold_entries = collect_entriesfor_physics(multiple_physics_toml, "cold")
  @test length(multiphys_heat_entries) == 1

  @test multiphys_cold_entries != multiphys_heat_entries


  multiple_arches_toml = TOML.parse("""
  [heat]
  cpu = ["test2"]
  cuda = ["cutest2", "final"]
  """)

  multiarches_entries = collect_entriesfor_physics(multiple_arches_toml, "heat")
  @test length(multiarches_entries) == 3
  @test allunique(multiarches_entries)
end

end
