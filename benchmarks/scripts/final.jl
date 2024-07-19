using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using DataFrames
using PrettyTables
using TOML
using JLD2

const slurm_id = ARGS[1]
const sim_name = ARGS[2]

# Collect meta config information to collect supported fields
# XXX: Add any extra headers here, config fields are automatically detected
# TODO: Have meta config information be in a seperate toml
meta_config_data = TOML.parsefile(get_config(sim_name, "cpu"))
meta_config = meta_config_data[string(0)]
const meta_field_names = split(meta_config["fields"], ",")
const solver_stages = split(meta_config["stages"], ",")

aggdatadir(args...) = datadir("exp_pro", sim_name, slurm_id, "autogen", args...)

function add_filename!(data_row, colname, filepath)
  filename = last(splitpath(filepath))
  push!(data_row, colname => filename)
end

function process_simdata(jld2file, jsonfile, config_data)
  data_row = Dict{String, Any}()

  benchmark_data = first(BenchmarkTools.load(jsonfile))
  task_key = only(keys(benchmark_data))

  push!(data_row, "Task ID" => task_key)
  add_filename!(data_row, "statsfile", jld2file)
  add_filename!(data_row, "benchfile", jsonfile)
  merge!(data_row, config_data[task_key]) # Adds sim parameters

  median_trial = median(benchmark_data[task_key])
  min_trial = minimum(benchmark_data[task_key])

  for stage in solver_stages
    push!(data_row, stage*" Median Time" => median_trial[stage].time)
    push!(data_row, stage*" Minimum Time" => min_trial[stage].time)

    push!(data_row, stage*" Median Mem" => median_trial[stage].memory)
    push!(data_row, stage*" Minimum Mem" => min_trial[stage].memory)
  end

  merge!(data_row, wload(jld2file)) # Adds result stats

  safesave(aggdatadir("results.jld2"), data_row)
end

function get_resultfiles(check)
  files = sort!(filter(filename -> occursin(check, filename), readdir(resultsdir(sim_name), join=true)))
end

function aggregate_data()
  for arch in ["cpu", "cuda"]
    jsonfiles = get_resultfiles("$(arch).json")
    jld2files = get_resultfiles("$(arch).jld2")

    config_name = get_config(sim_name, arch)
    isfile(config_name) || continue
    config_data = TOML.parsefile(config_name)

    for (jld2file, jsonfile) in zip(jld2files, jsonfiles)
      process_simdata(jld2file, jsonfile, config_data)
    end
  end
end

aggregate_data()

# TODO: Meant as a basic data processing pipeline
# Can create multiple scripts to roughly process data in general ways
pretty_results = collect_results(aggdatadir())

median_times = map(x -> x*" Median Time", solver_stages)
table_header = vcat(["Task ID", "statsfile", "benchfile"], meta_field_names, median_times, ["nf"])

select!(pretty_results, table_header)

for time in median_times
  transform!(pretty_results, [time] => ByRow(x -> x / 1e9) => [time])
end

# TODO: Can choose other sorting methods
for field_name in reverse(meta_field_names)
  sort!(pretty_results, Symbol(field_name))
end

mkpath(tablesdir(sim_name, slurm_id))

# TODO: Can change this backened to be different from markdown
open(tablesdir(sim_name, slurm_id, "final.md"), "w") do results_file
  conf = set_pt_conf(tf = tf_markdown)
  pretty_table_with_conf(conf, results_file, pretty_results; header = table_header)
end
