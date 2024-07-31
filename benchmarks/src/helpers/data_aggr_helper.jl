using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using DataFrames
using JLD2

function aggregate_data(slurm_id, sim_name)
  all_arches = get_supported_arches()
  for arch in all_arches
    config_name = get_config(sim_name, arch)
    isfile(config_name) || continue
    
    config_data = TOML.parsefile(config_name)
    num_entries = get_config_size(config_data)

    for task_id in 1:num_entries
      task_key = string(task_id)

      benchmark_filepath = get_benchfile(task_key, sim_name, arch)
      solve_stats_filepath = get_statsfile(task_key, sim_name, arch)
      if !isfile(benchmark_filepath) || !isfile(solve_stats_filepath)
        @info "Result file not found for $task_key on $arch, skipping"
        continue
      end

      aggregated_results = process_simdata(task_key, arch, benchmark_filepath, solve_stats_filepath, config_data)
      safesave(aggdatadir(sim_name, slurm_id, "results.jld2"), aggregated_results)
    end

  end
end

function process_simdata(task_key, arch, benchmark_filepath, solve_stats_filepath, config_data)
  data_row = Dict{String, Any}()

  add_debug_simdata!(data_row, task_key, arch)
  add_benchmark_data!(data_row, task_key, benchmark_filepath, config_data)
  add_solver_stats_data!(data_row, solve_stats_filepath)

  data_row
end

function add_debug_simdata!(data_row, task_key, arch)
  push!(data_row, "statsfile" => get_statsfile_name(task_key, arch))
  push!(data_row, "benchfile" => get_benchfile_name(task_key, arch))
  push!(data_row, "Task ID" => task_key)
end

function add_benchmark_data!(data_row, task_key, benchmark_filepath, config_data)
  benchmark_data = first(BenchmarkTools.load(benchmark_filepath))

  merge!(data_row, config_data[task_key]) # Adds sim parameters

  median_trial = median(benchmark_data[task_key])
  min_trial = minimum(benchmark_data[task_key])

  add_benchmark_data!(data_row, median_trial, "Median")
  add_benchmark_data!(data_row, min_trial, "Minimum")
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

function add_solver_stats_data!(data_row, solve_stats_filepath)
  merge!(data_row, wload(solve_stats_filepath))
end

