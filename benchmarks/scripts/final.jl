using DrWatson
@quickactivate "benchmarks"

using BenchmarkTools
using CSV
using DataFrames
using PrettyTables
using TOML
using JLD2

slurm_id = ARGS[1]
sim_name = ARGS[2]

simdir(args...) = srcdir(sim_name, args...)
paramsdir(args...) = simdir("params", args...)
resultsdir(args...) = datadir("sims", sim_name, args...)
tablesdir(args...) = datadir("exp_pro", sim_name, slurm_id, args...)

get_configfile_name(sim_name, code_target) = simdir("$(sim_name)_$(code_target).toml")
get_configfile(sim_name, code_target) = TOML.parsefile(get_configfile_name(sim_name, code_target))

# Collect meta config information to collect supported fields
# XXX: Add any extra headers here, config fields are automatically detected
# TODO: Have meta config information be in a seperate toml
meta_config_data = get_configfile(sim_name, "cpu")
meta_config = meta_config_data[string(0)]
meta_field_names = split(meta_config["fields"], ",")
solver_stages = split(meta_config["stages"], ",")

table_header = vcat(["Task ID"], meta_field_names, solver_stages, ["Steps", "Steps/Second"])
table_data = []

for code_target = ["cpu", "cuda"]

  # Collect all benchmark result files from "results" directory
  json_check = "$(code_target).json"
  jsonfiles = sort!(filter(filename -> occursin(json_check, filename), readdir(resultsdir(), join=true)))

  jld2_check = "$(code_target).jld2"
  jld2files = sort!(filter(filename -> occursin(jld2_check, filename), readdir(resultsdir(), join=true)))

  # Collect config information for indvidual task configs
  config_name = get_configfile_name(sim_name, code_target)

  if !isfile(config_name)
    continue
  end

  config_data = get_configfile(sim_name, code_target)

  # Read through all solution files and collect data into Matrix
  for (jld2file, jsonfile) in zip(jld2files, jsonfiles)
    data_row = []

    # Get ID of task from benchmark group name
    benchmark_data = first(BenchmarkTools.load(jsonfile))
    task_key = only(keys(benchmark_data))
    push!(data_row, task_key)

    # Collect config information to display in table
    curr_config = config_data[task_key]
    for field_name in meta_field_names
      push!(data_row, curr_config[field_name])
    end

    # TODO: Median can be swapped out for any other statistics
    trial = median(benchmark_data[task_key])

    # TODO: Can extend this to also collect memory
    for stage in solver_stages
      push!(data_row, trial[stage].time / 1e9)
    end

    result_stats_data = wload(jld2file)
    steps_taken = result_stats_data["naccept"]
    push!(data_row, steps_taken)

    # TODO: This still checks for "Solve" stage
    step_per_second = floor(Int64, steps_taken / (trial["Solve"].time / 1e9))
    push!(data_row, step_per_second)

    push!(table_data, data_row)
  end
end

if length(table_data) <= 1
  error("Data output is too short, please include more than one result")
end

# TODO: Need to permute the Matrix in order to add in data to DataFrame as cols
table = permutedims(foldl(hcat, table_data), (2, 1))

num_cols = length(table_header)

# Add data to DataFrame
data_frame = DataFrame()
for (i, colname) in enumerate(table_header)
  data_frame[!, colname] = table[:, i]
end

# Sort by all config information, with the first fields of a table having higher priority
# TODO: Can choose other sorting methods
for field_name in reverse(meta_field_names)
  sort!(data_frame, Symbol(field_name))
end

mkpath(tablesdir())
# Can be used to save table for later use
CSV.write(tablesdir("final.csv"), data_frame)
# csv_table = CSV.read("final.csv", DataFrame)

# TODO: Can change this backened to be different from markdown
# Open file where table will be saved
open(tablesdir("final.md"), "w") do results_file
  conf = set_pt_conf(tf = tf_markdown)
  pretty_table_with_conf(conf, results_file, data_frame; header = table_header)
end
