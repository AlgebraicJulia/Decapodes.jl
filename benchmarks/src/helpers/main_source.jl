using DrWatson
@quickactivate :benchmarks

using TOML
using MLStyle

include(helpersdir("physics_config_helper.jl"))

const DEFAULT_CPU_SLURMARGS = `--partition=hpg-milan`
const DEFAULT_CUDA_SLURMARGS = `--partition=gpu --gres=gpu:a100:1`
const DEFAULT_MAXJOBS = 6

function generate_all_configs()
  for physics in listof_main_physics()
      generate_configsfor_physics(physics)
  end
end

function run_all_physics()
  main_config_info = load_main_config()
  validate_all_physics(main_config_info)

  for physics in listof_main_physics()
      physics_configs = collect_simsfor_physics(main_config_info, physics)
      run_single_physics(physics, physics_configs)
  end
end

function validate_all_physics(main_config_info)
  for physics in listof_main_physics()
      physics_configs = collect_simsfor_physics(main_config_info, physics)
      for config in physics_configs
          is_valid_config_instance(config)
      end
  end
end

function is_valid_config_instance(sim_namedata::SimNameData)
  is_supported_arch(sim_namedata.arch) || error("Architecture $(arch) is not in list $(supported_arches())")
  simfile = physicsfile_path(sim_namedata)
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

  !isempty(dependency_ids) || return

  run(`sbatch --dependency=afterok:$(join(dependency_ids, ",")) $(scriptsdir("final.sh")) $physics`)
end

function run_single_config!(dependency_ids, sim_namedata::SimNameData)

  main_config_info = load_main_config()
  concurrent_jobs = get_concur_jobs(main_config_info, sim_namedata)
  slurm_args = get_slurm_args(main_config_info, sim_namedata)

  physics = sim_namedata.physics
  arch = sim_namedata.arch
  tag = sim_namedata.tag

  count = simconfig_size(load_simconfig(sim_namedata))

  @info "Running $count $arch tasks"

  jobid = readchomp(`sbatch --array=1-$(count)%$(concurrent_jobs) --parsable $(slurm_args) $(scriptsdir("array.sh")) $physics $arch $tag`)

  @info "Job ID for $(tagfor_run(sim_namedata)) is: $jobid"
  push!(dependency_ids, jobid)
  return dependency_ids
end

function get_slurm_args(main_config_info, snd::SimNameData)

  slurm_args = get_config_arg(main_config_info, snd, "slurm_args")

  if isnothing(slurm_args)
    return @match snd.arch begin
      "cpu" => DEFAULT_CPU_SLURMARGS
      "cuda" => DEFAULT_CUDA_SLURMARGS
    end
  end

  return Cmd(slurm_args)
end

function get_concur_jobs(main_config_info, snd::SimNameData)
  job_arg = get_config_arg(main_config_info, snd, "concur_jobs")
  return isnothing(job_arg) ? DEFAULT_MAXJOBS : job_arg
end
