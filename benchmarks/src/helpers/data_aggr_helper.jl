using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using DataFrames
using JLD2

function aggregate_data(slurm_id, physics)

  main_config_info = load_main_config()
  sims_to_process = collect_mainconfig_simentries(physics, main_config_info)
  for sim_namedata in sims_to_process

    # Meant to work when not all sims in main_config are run
    if !all_resultfiles_exist(sim_namedata)
      @info "Result file not found for $(sim_namedata), skipping"
      continue
    end

    benchmark_filepath = benchfile_path(sim_namedata)
    solve_stats_filepath = statsfile_path(sim_namedata)
    config_data = TOML.parsefile(simconfig_path(sim_namedata))

    aggregated_results = process_simdata(sim_namedata, benchmark_filepath, solve_stats_filepath, config_data)
    safesave(aggdatadir(physics, slurm_id, "results.jld2"), aggregated_results)

  end
end

function all_resultfiles_exist(sim_namedata)
  return isfile(benchfile_path(sim_namedata)) && isfile(statsfile_path(sim_namedata))
end

# TODO: Come back and clean up this function
function collect_mainconfig_simentries(physics, main_config_info)

  entries = SimNameData[]

  physics_configurations = collect_simsfor_physics(main_config_info, physics)
  for tagged_sim_namedata in physics_configurations

      tagged_sim_config = simconfig_path(tagged_sim_namedata)
      if !isfile(tagged_sim_config)
          @info "Config file for $tagged_sim_namedata not found, skipping"
          continue
      end

      tagged_sim_data = TOML.parsefile(tagged_sim_config)
      num_entries = autoconfig_size(tagged_sim_data)
      for task_id in 1:num_entries

          physics = tagged_sim_namedata.physics
          arch = tagged_sim_namedata.arch
          tag = tagged_sim_namedata.tag
          task_key = string(task_id)

          push!(entries, SimNameData(physics, arch, tag, task_key))
      end
  end

  entries
end

function process_simdata(sim_namedata::SimNameData, benchmark_filepath, solve_stats_filepath, config_data)
  data_row = Dict{String, Any}()

  add_debug_simdata!(data_row, sim_namedata)
  add_config_data!(data_row, sim_namedata.task_key, config_data)

  benchmark_data = first(BenchmarkTools.load(benchmark_filepath))
  add_benchmark_data!(data_row, sim_namedata.task_key, benchmark_data)

  stats_data = wload(solve_stats_filepath)
  add_solver_stats_data!(data_row, stats_data)

  data_row
end

function add_debug_simdata!(data_row, sim_namedata::SimNameData)
  push!(data_row, "statsfile" => statsfile_name(sim_namedata))
  push!(data_row, "benchfile" => benchfile_name(sim_namedata))
  push!(data_row, "Task ID" => sim_namedata.task_key)
  push!(data_row, "Tagged Name" => sim_namedata.tag)
  push!(data_row, "Architecture" => sim_namedata.arch)
  push!(data_row, "Simulation Name" => sim_namedata.physics)
end

function add_config_data!(data_row, task_key, config_data)
  merge!(data_row, config_data[task_key]) # Adds sim parameters
end

function add_benchmark_data!(data_row, task_key, benchmark_data)
  median_trial = median(benchmark_data[task_key])
  min_trial = minimum(benchmark_data[task_key])

  add_trial_data!(data_row, median_trial, "Median")
  add_trial_data!(data_row, min_trial, "Minimum")
end

function add_trial_data!(data_row, trial_data, name::String)
  stages = solver_stages()
  for stage in stages
    push!(data_row, get_benchmark_headername(stage, name, "Time") => trial_data[stage].time)
    push!(data_row, get_benchmark_headername(stage, name, "Mem") => trial_data[stage].memory)
  end
end

function get_benchmark_headername(stage::String, name::String, category::String)
  return "$(stage) $(name) $(category)"
end

function add_solver_stats_data!(data_row, stats_data)
  merge!(data_row, stats_data)
end
