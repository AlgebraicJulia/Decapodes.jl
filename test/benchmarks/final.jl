using BenchmarkTools
using TOML
using PrettyTables
using DataFrames

println("Combining all files")

json_check = r"\.json$"
jsonfiles = filter(filename -> occursin(json_check, filename), readdir("results", join=true))

all_config_data = TOML.parsefile("benchmarks_config_heat.toml")

meta_config = all_config_data[string(0)]
meta_field_names = split(meta_config["fields"], ",")

results_file = open("final.md", "w")

table_header = vcat(meta_field_names, ["Setup(s)", "Mesh(s)", "Simulate(s)", "Solve(s)"])
table_data = []

for jsonfile in jsonfiles
  data_row = []
  benchmark_data = first(BenchmarkTools.load(jsonfile))
  task_key = only(keys(benchmark_data))
  trial = median(benchmark_data[task_key])

  curr_config = all_config_data[task_key]
  for field_name in meta_field_names
    push!(data_row, curr_config[field_name])
  end
  for stage in keys(trial)
    push!(data_row, trial[stage].time / 1e9)
  end
  push!(table_data, data_row)
end

table = permutedims(foldl(hcat, table_data), (2, 1))

num_cols = length(table_header)

data_frame = DataFrame()
for (i, colname) in enumerate(table_header)
  data_frame[!, colname] = table[:, i]
end

for field_name in reverse(meta_field_names)
  sort!(data_frame, Symbol(field_name))
end

conf = set_pt_conf(tf = tf_markdown)
pretty_table_with_conf(conf, results_file, data_frame; header = table_header)

close(results_file)
