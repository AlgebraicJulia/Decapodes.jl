# Collection of path function mimicking DrWatson's functions
using DrWatson
@quickactivate :benchmarks

export resultsdir, tablesdir, helpersdir, aggdatadir, postprocessdir
export get_configname, get_config, get_simfile_name, get_simfile, get_statsfile_name, get_statsfile, get_benchfile_name, get_benchfile
export get_config_size, get_meta_config_info

helpersdir(args...) = srcdir("helpers", args...)

resultsdir(sim_name, args...) = datadir("sims", sim_name, args...)

tablesdir(sim_name, slurm_id, args...) = datadir("exp_pro", sim_name, slurm_id, args...)
aggdatadir(sim_name, slurm_id, args...) = tablesdir(sim_name, slurm_id, "autogen", args...)
postprocessdir(args...) = scriptsdir("post_processing", args...)

function get_configname(sim_name, arch)
  return "$(sim_name)_$(arch).toml"
end

function get_config(sim_name, arch)
  return srcdir(sim_name, get_configname(sim_name, arch))
end

get_simfile_name(sim_name) = return "$sim_name.jl"

function get_simfile(sim_name)
  return srcdir(sim_name, get_simfile_name(sim_name))
end

get_statsfile_name(task_key, arch) = return "stats_$(task_key)_$(arch).jld2"

function get_statsfile(task_key, sim_name, arch)
  return resultsdir(sim_name, get_statsfile_name(task_key, arch))
end

get_benchfile_name(task_key, arch) = return "benchmarks_$(task_key)_$(arch).json"

function get_benchfile(task_key, sim_name, arch)
  return resultsdir(sim_name, get_benchfile_name(task_key, arch))
end

function get_config_size(config_data)
  key_list = keys(config_data)
  return length(key_list) - 1 # Don't include meta info
end

function get_meta_config_info(benchmark_config)
  return benchmark_config[meta_config_id()]
end