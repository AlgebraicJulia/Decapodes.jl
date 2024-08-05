using DrWatson
@quickactivate :benchmarks

include(helpersdir("config_helper.jl"))

using TOML
using MLStyle

const default_config_file = srcdir("config.toml")

if length(ARGS) == 0
  config_file = default_config_file
elseif length(ARGS) == 1
  config_file = ARGS[1]
else
  error("Usage: [config file name]")
end

process_benchmark_config(config_file)


