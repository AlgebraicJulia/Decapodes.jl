using DrWatson
@quickactivate "benchmarks"

export physicsdir, resultsdir, tablesdir, helpersdir, aggdatadir, postprocessdir
export tagfor_run, tagfor_task
export simconfig_name, simconfig_path, physics_filename, physicspath, physicsconfig_path, statsfile_name, statsfile_path, benchfile_name, benchfile_path
export autoconfig_size, get_meta_config_info

export SimNameData
import Base.show

helpersdir(args...) = srcdir("helpers", args...)

physicsdir(args...) = srcdir("physics", args...)

resultsdir(physics, args...) = datadir("sims", physics, args...)

tablesdir(physics, slurm_id, args...) = datadir("exp_pro", physics, slurm_id, args...)
aggdatadir(physics, slurm_id, args...) = tablesdir(physics, slurm_id, "autogen", args...)
postprocessdir(args...) = scriptsdir("post_processing", args...)

struct SimNameData
  physics::String
  arch::String
  tag::String
  task_key::String
end

function SimNameData(physics::String, arch::String, tag::String)
  return SimNameData(physics, arch, tag, "0")
end

function SimNameData(physics::String, arch::String, tag::String, task_id::Int)
  return SimNameData(physics, arch, tag, string(task_id))
end

function Base.show(io::IO, snd::SimNameData)
  if snd.task_key == "0"
    print(io, "$(snd.physics) on $(snd.arch) tagged as $(snd.tag)")
  else
    print(io, "$(snd.physics) on $(snd.arch) tagged as $(snd.tag) running for task $(snd.task_key)")
  end
end

function tagfor_run(simdata::SimNameData)
  return "$(simdata.physics)_$(simdata.arch)_$(simdata.tag)"
end

function tagfor_task(simdata)
  return "$(tagfor_run(simdata))_$(simdata.task_key)"
end

function simconfig_name(simdata::SimNameData)
  return "$(tagfor_run(simdata)).toml"
end

function simconfig_path(simdata::SimNameData)
  return physicsdir(simdata.physics, simconfig_name(simdata))
end

physics_filename(physics) = return "$physics.jl"

function physicspath(physics)
  return physicsdir(physics, physics_filename(physics))
end

function physicsconfig_path(physics)
  physicsdir(physics, physics_config_filename())
end

statsfile_name(simdata::SimNameData) = return "stats_$(tagfor_task(simdata)).jld2"

function statsfile_path(simdata::SimNameData)
  return resultsdir(simdata.physics, statsfile_name(simdata))
end

benchfile_name(simdata::SimNameData) = return "benchmarks_$(tagfor_task(simdata)).json"

function benchfile_path(simdata::SimNameData)
  return resultsdir(simdata.physics, benchfile_name(simdata))
end

function autoconfig_size(config_data)
  key_list = keys(config_data)
  return length(key_list) - 1 # Don't include meta info
end

function get_meta_config_info(benchmark_config)
  return benchmark_config[meta_config_id()]
end
