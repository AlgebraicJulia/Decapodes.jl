# Collection of path function mimicking DrWatson's functions
using DrWatson

export paramsdir, resultsdir, tablesdir, get_configname, 
    get_config, get_simfile, get_statsfile, get_benchfile

paramsdir(sim_name, args...) = srcdir(sim_name, "params", args...)
resultsdir(sim_name, args...) = datadir("sims", sim_name, args...)
tablesdir(sim_name, slurm_id, args...) = datadir("exp_pro", sim_name, slurm_id, args...)

function get_configname(sim_name, arch)
    return "$(sim_name)_$(arch).toml"
end

function get_config(sim_name, arch)
    return srcdir(sim_name, get_configname(sim_name, arch))
end

function get_simfile(sim_name)
    return srcdir(sim_name, "$sim_name.jl")
end

function get_statsfile(task_key, sim_name, arch)
    return resultsdir(sim_name, "stats_$(task_key)_$(arch).jld2")
end

function get_benchfile(task_key, sim_name, arch)
    return resultsdir(sim_name, "benchmarks_$(task_key)_$(arch).json")
end