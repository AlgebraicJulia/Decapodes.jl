using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using DataFrames
using JLD2

function aggregate_data(slurm_id, sim_name)

  main_config_info = TOML.parsefile(mainconfig_path())
  sims_to_process = collect_mainconfig_simentries(sim_name, main_config_info)
  for sim_namedata in sims_to_process

    benchmark_filepath = benchfile_path(sim_namedata)
    solve_stats_filepath = statsfile_path(sim_namedata)
    if !isfile(benchmark_filepath) || !isfile(solve_stats_filepath)
      @info "Result file not found for $(sim_namedata.task_key) on $(sim_namedata.arch), skipping"
      continue
    end

    config_data = TOML.parsefile(config_path(sim_namedata))

    aggregated_results = process_simdata(sim_namedata, benchmark_filepath, solve_stats_filepath, config_data)
    safesave(aggdatadir(sim_name, slurm_id, "results.jld2"), aggregated_results)

  end
end

function collect_mainconfig_simentries(sim_name, main_config_info)
  
  haskey(main_config_info, sim_name) || error("Main config missing $(sim_name) entry")
  sim_config_info = main_config_info[sim_name]

  entries = SimNameData[]
  for arch in keys(sim_config_info)
      for tag in sim_config_info[arch]

          tagged_sim_namedata = SimNameData(sim_name, arch, tag)
          tagged_sim_config = config_path(tagged_sim_namedata)
          if !isfile(tagged_sim_config) 
            @info "Config file for $sim_name on $arch named $tag not found, skipping"
            continue
          end
      
          tagged_sim_data = TOML.parsefile(tagged_sim_config)
          num_entries = get_autoconfig_size(tagged_sim_data)
          for task_id in 1:num_entries
              task_key = string(task_id)
              push!(entries, SimNameData(sim_name, arch, tag, task_key))
          end
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
  push!(data_row, "Simulation Name" => sim_namedata.sim_name)
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

