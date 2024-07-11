using BenchmarkTools
using TOML
using CSV
using PrettyTables
using DataFrames

save_dir = ARGS[1]
sim_name = ARGS[2]

cd(sim_name)

# Collect all benchmark result files from "results" directory
json_check = r"\.json$"
jsonfiles = sort!(filter(filename -> occursin(json_check, filename), readdir("results", join=true)))

txt_check = r"\.txt$"
txtfiles = sort!(filter(filename -> occursin(txt_check, filename), readdir("results", join=true)))

# Collect config information for indvidual task configs
all_config_data = TOML.parsefile("$(sim_name)_cpu.toml")

# Collect meta config information to collect supported fields
meta_config = all_config_data[string(0)]
meta_field_names = split(meta_config["fields"], ",")
solver_stages = split(meta_config["stages"], ",")

# XXX: Add any extra headers here, config fields are automatically detected
table_header = vcat(["Task ID"], meta_field_names, solver_stages, ["Steps", "Steps/Second"])
table_data = []

# Read through all solution files and collect data into Matrix
for (txtfile, jsonfile) in zip(txtfiles, jsonfiles)
  data_row = []

  # Get ID of task from benchmark group name
  benchmark_data = first(BenchmarkTools.load(jsonfile))
  task_key = only(keys(benchmark_data))
  push!(data_row, task_key)

  # Collect config information to display in table
  curr_config = all_config_data[task_key]
  for field_name in meta_field_names
    push!(data_row, curr_config[field_name])
  end

  # TODO: Median can be swapped out for any other statistics
  trial = median(benchmark_data[task_key])

  # TODO: Can extend this to also collect memory
  for stage in solver_stages
    push!(data_row, trial[stage].time / 1e9)
  end

  stats_data = split(readline(txtfile), ",")
  steps_taken = parse(Int64, stats_data[11]) # Assume DEStats is constant here
  push!(data_row, steps_taken)

  # TODO: This still checks for "Solve" stage
  step_per_second = floor(Int64, steps_taken / (trial["Solve"].time / 1e9))
  push!(data_row, step_per_second)

  push!(table_data, data_row)
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

mkpath(save_dir)
cd(save_dir)

# Can be used to save table for later use
CSV.write("final.csv", data_frame)
# csv_table = CSV.read("final.csv", DataFrame)

# TODO: Can change this backened to be different from markdown
# Open file where table will be saved
open("final.md", "w") do results_file
  conf = set_pt_conf(tf = tf_markdown)
  pretty_table_with_conf(conf, results_file, data_frame; header = table_header)
end
