# Collection of path function mimicking DrWatson's functions
using DrWatson
@quickactivate :benchmarks

export physicsdir, resultsdir, tablesdir, helpersdir, aggdatadir, postprocessdir
export SimNameData, tagfor_run, tagfor_task
export config_name, config_path, get_simfile_name, get_simfile, statsfile_name, statsfile_path, benchfile_name, benchfile_path, mainconfig_path
export get_autoconfig_size, get_meta_config_info

helpersdir(args...) = srcdir("helpers", args...)

physicsdir(args...) = srcdir("physics", args...)

resultsdir(sim_name, args...) = datadir("sims", sim_name, args...)

tablesdir(sim_name, slurm_id, args...) = datadir("exp_pro", sim_name, slurm_id, args...)
aggdatadir(sim_name, slurm_id, args...) = tablesdir(sim_name, slurm_id, "autogen", args...)
postprocessdir(args...) = scriptsdir("post_processing", args...)

struct SimNameData
  sim_name::String
  arch::String
  tag::String
  task_key::String
end

function SimNameData(sim_name::String, arch::String, tag::String)
  return SimNameData(sim_name, arch, tag, "0")
end


function tagfor_run(simdata::SimNameData)
  return "$(simdata.sim_name)_$(simdata.arch)_$(simdata.tag)"
end

function tagfor_task(simdata)
  return "$(tagfor_run(simdata))_$(simdata.task_key)"
end


function config_name(simdata::SimNameData)
  return "$(tagfor_run(simdata)).toml"
end

function config_path(simdata::SimNameData)
  return physicsdir(simdata.sim_name, config_name(simdata))
end

get_simfile_name(sim_name) = return "$sim_name.jl"

function get_simfile(sim_name)
  return physicsdir(sim_name, get_simfile_name(sim_name))
end

mainconfig_path() = srcdir("main_config.toml")

statsfile_name(simdata::SimNameData) = return "stats_$(tagfor_task(simdata)).jld2"

function statsfile_path(simdata::SimNameData)
  return resultsdir(simdata.sim_name, statsfile_name(simdata))
end

benchfile_name(simdata::SimNameData) = return "benchmarks_$(tagfor_task(simdata)).json"

function benchfile_path(simdata::SimNameData)
  return resultsdir(simdata.sim_name, benchfile_name(simdata))
end

function get_autoconfig_size(config_data)
  key_list = keys(config_data)
  return length(key_list) - 1 # Don't include meta info
end

function get_meta_config_info(benchmark_config)
  return benchmark_config[meta_config_id()]
end