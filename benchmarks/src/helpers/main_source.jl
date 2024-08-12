using DrWatson
@quickactivate :benchmarks

using TOML
using MLStyle

include(helpersdir("main_config_helper.jl"))
include(helpersdir("config_helper.jl"))

function generate_all_configs()
  for physics in listof_main_physics()
      process_benchmark_config(physicsconfig_path(physics))
  end
end

function run_all_physics()
  main_config_info = TOML.parsefile(mainsim_config_path())
  validate_all_physics(main_config_info)

  for physics in listof_main_physics()
      physics_configs = collect_entriesfor_physics(main_config_info, physics)
      run_single_physics(physics, physics_configs)
  end
end

function listof_main_physics()
  main_config_info = TOML.parsefile(mainsim_config_path())
  return collect(keys(main_config_info))
end

function validate_all_physics(main_config_info)
  for physics in listof_main_physics()
      physics_configs = collect_entriesfor_physics(main_config_info, physics)
      for config in physics_configs
          validate_config_instance(config)
      end
  end
end

function validate_config_instance(sim_namedata::SimNameData)
  physics = sim_namedata.physics

  is_supported_arch(sim_namedata.arch) || error("Architecture $(arch) is not in list $(supported_arches())")
  simfile = physicspath(physics)
  if !isfile(simfile)
      error("Simulation file at $(simfile) was not found")
  end

  config = simconfig_path(sim_namedata)
  if !isfile(config)
      error("Config file at $(config) was not found")
  end
  return true
end

function run_single_physics(physics, physics_configs)

  if isempty(physics_configs)
      return
  end

  rm(resultsdir(physics), recursive=true, force=true)
  mkpath(resultsdir(physics))

  dependency_ids = []

  for config in physics_configs
      run_single_config!(dependency_ids, config)
  end

  if isempty(dependency_ids)
      return
  end

  run(`sbatch --dependency=afterok:$(join(dependency_ids, ",")) $(scriptsdir("final.sh")) $physics`)
end

function run_single_config!(dependency_ids, sim_namedata::SimNameData)
  config_data = TOML.parsefile(simconfig_path(sim_namedata))

  concurrent_jobs = 6

  physics = sim_namedata.physics
  arch = sim_namedata.arch
  tag = sim_namedata.tag

  count = autoconfig_size(config_data)

  @info "Running $count $arch tasks"

  args = @match arch begin
    "cpu" => `--partition=hpg-milan`
    "cuda" => `--partition=gpu --gres=gpu:a100:1`
  end

  jobid = readchomp(`sbatch --array=1-$(count)%$(concurrent_jobs) --parsable $(args) $(scriptsdir("array.sh")) $physics $arch $tag`)

  @info "Job ID for $(tagfor_run(sim_namedata)) is: $jobid"
  push!(dependency_ids, jobid)
end
