module BenchmarkConfig

using TOML
using CSV
using MLStyle
using Base.Iterators

const metaconfig_file = "metaconfig.csv"
const dirty_bit = "dirty.txt"

function shorten_code_target(target)
  @match target begin
    "CPUTarget" => "cpu"
    "CUDATarget" => "cuda"
    _ => error("Please supply either 'CPUTarget' or 'CUDATarget'")
  end
end

function parse_float_type(float_type)
  @match float_data begin
    "Float32" => Float32
    "Float64" => Float64
    "ComplexF32" => ComplexF32
    "ComplexF64" => ComplexF64
    _ => error("Float data $(float_data) is not valid, exiting early") 
  end
end

function validate_metaconfig(metaconfig_file)
  if !isfile(metaconfig_file)
    error("Please create a file '$metaconfig_file'")
  end

  csv_config = CSV.File(metaconfig_file)

  temp_key_store = csv_config.names

  if length(temp_key_store) < 2 || temp_key_store[1] != :sim_name || temp_key_store[2] != :code_target
    error("Improper field names, please begin with 'sim_name' and 'code_target'")
  end

  if isempty(csv_config)
    error("'$metaconfig_file' is empty, please populate it")
  end
end

function generate_config(csv_row)

  target_name = shorten_code_target(csv_row.code_target)

  # TODO: Move this further out
  stages = ["Setup", "Mesh", "Simulate", "Solve"]

  config_store = []
  for csv_val in csv_row
    temp_split = split(csv_val)
    if isempty(temp_split)
      error("Config file found an empty value in $(metaconfig_file), exiting...")
    end
    push!(config_store, temp_split)
  end

  data = Dict()

  for (idx, config) in enumerate(product(config_store...))
      datum = Dict()
      for (col_idx, val) in enumerate(config)
        # Offset to avoid reading sim_name
        push!(datum, keys(csv_row)[col_idx] => val)
      end
      push!(data, string(idx) => datum)
  end

  open("count_$(target_name).txt", "w") do file
    write(file, string(length(data)))
  end

    # Avoid reading in sim name
  meta_data = Dict("fields" => join(keys(csv_row)[2:end], ","),
    "stages" => join(stages, ","))
  push!(data, string(0) => meta_data)

  open("$(csv_row.sim_name)_$(target_name).toml", "w") do io
    TOML.print(io, data)
  end
end

function update_tracker(file_name)
  open(file_name, "w") do file
    write(file, "Configuration changed since last run, 
                please reset stored benchmarking parameters, then delete me")
  end
end

function run_metaconfig(metaconfig_file)
  validate_metaconfig(metaconfig_file)

  csv_config = CSV.File(metaconfig_file)

  for csv_row in csv_config
    cd(csv_row.sim_name) do
      generate_config(csv_row)
      update_tracker(dirty_bit)
    end
  end
end

run_metaconfig(metaconfig_file)

end