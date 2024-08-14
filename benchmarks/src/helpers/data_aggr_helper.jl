using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using DataFrames
using JLD2

function aggregate_data(slurm_id, physics)
  sims_to_process = collect_mainconfig_simentries(physics)

  for sim_namedata in sims_to_process

    # Meant to work when not all sims in main_config are run
    if !all_resultfiles_exist(sim_namedata)
      @info "Result file not found for $(sim_namedata), skipping"
      continue
    end

    aggregated_results = process_simdata(sim_namedata)

    savefile_path = aggdatadir(physics, slurm_id, "results.jld2")
    safesave(savefile_path, aggregated_results)

  end
end

function all_resultfiles_exist(sim_namedata)
  return isfile(benchfile_path(sim_namedata)) && isfile(statsfile_path(sim_namedata))
end

# TODO: Come back and clean up this function
function collect_mainconfig_simentries(physics)

  entries = SimNameData[]

  physics_configurations = collect_simsfor_physics(physics)
  for sim_namedata in physics_configurations

      if !isfile(simconfig_path(sim_namedata))
          @info "Config file for $sim_namedata not found, skipping"
          continue
      end

      config_entries = simconfig_size(load_simconfig(sim_namedata))
      for task_id in 1:config_entries

          physics = sim_namedata.physics
          arch = sim_namedata.arch
          tag = sim_namedata.tag
          task_key = string(task_id)

          push!(entries, SimNameData(physics, arch, tag, task_key))
      end
  end

  entries
end

function process_simdata(snd::SimNameData)
  return process_simdata(snd, load_benchfile(snd), load_statsfile(snd), load_simconfig(snd))
end

function process_simdata(sim_namedata::SimNameData, benchmark_data, stats_data, config_data)
  data_row = Dict{String, Any}()

  add_debug_simdata!(data_row, sim_namedata)
  add_config_data!(data_row, sim_namedata.task_key, config_data)
  add_benchmark_data!(data_row, sim_namedata.task_key, benchmark_data)
  add_solver_stats_data!(data_row, stats_data)

  return data_row
end

function add_debug_simdata!(data_row, sim_namedata::SimNameData)
  push!(data_row, "statsfile" => statsfile_name(sim_namedata))
  push!(data_row, "benchfile" => benchfile_name(sim_namedata))
  push!(data_row, "Task ID" => sim_namedata.task_key)
  push!(data_row, "Tagged Name" => sim_namedata.tag)
  push!(data_row, "Architecture" => sim_namedata.arch)
  push!(data_row, "Simulation Name" => sim_namedata.physics)
end

function add_config_data!(data_row, task_key::String, config_data)
  merge!(data_row, config_data[task_key]) # Adds sim parameters
end

function add_benchmark_data!(data_row, task_key::String, benchmark_data)
  median_trial = median(benchmark_data[task_key])
  min_trial = minimum(benchmark_data[task_key])

  add_trial_data!(data_row, "Median", median_trial)
  add_trial_data!(data_row, "Minimum", min_trial)
end

# TODO: Add benchmark parameters for each stage
function add_trial_data!(data_row, statistic_name::String, trial_data)
  for stage in solver_stages()
    stage_data = trial_data[stage]
    for field in fieldnames(typeof(stage_data))
      field != :params || continue
      push!(data_row, get_benchmark_headername(stage, statistic_name, String(field)) => getfield(stage_data, field))
    end
  end
end

function get_benchmark_headername(stage::String, statistic_name::String, category::String)
  return "$(stage) $(statistic_name) $(category)"
end

function add_solver_stats_data!(data_row, stats_data)
  merge!(data_row, stats_data)
end
