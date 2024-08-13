using DrWatson
@quickactivate "benchmarks"

using TOML

export physicsdir, resultsdir, tablesdir, helpersdir, aggdatadir, postprocessdir
export tagfor_run, tagfor_task
export simconfig_name, simconfig_path, physicsfile_name, physicsfile_path, physicsconfig_path, physicsconfig_name, statsfile_name, statsfile_path, benchfile_name, benchfile_path
export load_simconfig, load_physicsconfig
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
    print(io, "$(snd.physics) on $(snd.arch) tagged as '$(snd.tag)'")
  else
    print(io, "$(snd.physics) on $(snd.arch) tagged as '$(snd.tag)' running for task $(snd.task_key)")
  end
end

function tagfor_run(simdata::SimNameData)
  return "$(simdata.physics)_$(simdata.arch)_$(simdata.tag)"
end

function tagfor_task(simdata::SimNameData)
  return "$(tagfor_run(simdata))_$(simdata.task_key)"
end

function simconfig_name(simdata::SimNameData)
  return "$(tagfor_run(simdata)).toml"
end

load_simconfig(simdata::SimNameData) = TOML.parsefile(simconfig_path(simdata))
simconfig_path(simdata::SimNameData) = physicsdir(simdata.physics, simconfig_name(simdata))

physicsfile_name(simdata::SimNameData) = return "$(simdata.physics).jl"
physicsfile_path(simdata::SimNameData) = physicsdir(simdata.physics, physicsfile_name(simdata))

load_physicsconfig(simdata::SimNameData) = TOML.parsefile(physicsconfig_path(simdata))
physicsconfig_name() = return "config.toml"
physicsconfig_path(simdata::SimNameData) = physicsdir(simdata.physics, physicsconfig_name())


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
