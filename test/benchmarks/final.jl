using BenchmarkTools
using TOML
println("Combining all files")

json_check = r"\.json$"
jsonfiles = filter(filename -> occursin(json_check, filename), readdir("results", join=true))

all_config_data = TOML.parsefile("benchmarks_config_heat.toml")

results_file = open("final.log", "w")

for jsonfile in jsonfiles
  benchmark_data = BenchmarkTools.load(jsonfile)[begin]
  task_key = only(keys(benchmark_data))
  trial = median(benchmark_data[task_key])
  write(results_file, string(all_config_data[task_key])*"\n")
  for stage in keys(trial)
    write(results_file, "$(stage)| Time(s): $(trial[stage].time / 1e9), Memory(GB): $(trial[stage].memory / 1e9)\n")
  end
  write(results_file, "\n")
end

close(results_file)
