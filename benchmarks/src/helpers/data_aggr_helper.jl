using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using DataFrames
using JLD2

function aggregate_data(slurm_id, sim_name)
  for arch in ["cpu", "cuda"]
    config_name = get_config(sim_name, arch)
    isfile(config_name) || continue
    
    config_data = TOML.parsefile(config_name)

    jsonfiles = get_resultfiles("$(arch).json")
    jld2files = get_resultfiles("$(arch).jld2")

    for (jld2file, jsonfile) in zip(jld2files, jsonfiles)
      data_row = process_simdata(jld2file, jsonfile, config_data)
      safesave(aggdatadir(sim_name, slurm_id, "results.jld2"), data_row)
    end
  end
end

function get_resultfiles(filename_check)
  resultsdir_filepaths = readdir(resultsdir(sim_name), join=true)
  sort!(filter(filename -> occursin(filename_check, filename), resultsdir_filepaths))
end

function process_simdata(jld2file, jsonfile, config_data)
  data_row = Dict{String, Any}()

  # Passed in file name for debug info
  add_filename!(data_row, jld2file, "statsfile")
  add_filename!(data_row, jsonfile, "benchfile")

  # Adds data from benchmark
  benchmark_data = first(BenchmarkTools.load(jsonfile))
  task_key = only(keys(benchmark_data))

  push!(data_row, "Task ID" => task_key)
  merge!(data_row, config_data[task_key]) # Adds sim parameters

  median_trial = median(benchmark_data[task_key])
  min_trial = minimum(benchmark_data[task_key])

  add_benchmark_data!(data_row, median_trial, "Median")
  add_benchmark_data!(data_row, min_trial, "Minimum")

  # Adds data from solver statistics
  merge!(data_row, wload(jld2file))

  data_row
end

function add_filename!(data_row, filepath, colname::String)
  filename = last(splitpath(filepath))
  push!(data_row, colname => filename)
end

function add_benchmark_data!(data_row, trial_data, name::String)
  solver_stages = get_solver_stages()
  for stage in solver_stages
    push!(data_row, get_benchmark_headername(stage, name, "Time") => trial_data[stage].time)
    push!(data_row, get_benchmark_headername(stage, name, "Mem") => trial_data[stage].memory)
  end
end

function get_benchmark_headername(stage::String, name::String, category::String)
  return stage*" $(name) $(category)"
end