using DrWatson
@quickactivate :benchmarks
using TOML

function process_benchmark_config(configname)
    benchmark_config = TOML.parsefile(configname)
    validate_config(benchmark_config)
    load_save_benchmark_data(benchmark_config)
end

# Keep validation seperate to prevent creation of intermediate files/folders if error
function validate_config(benchmark_config)
  if isempty(benchmark_config)
    error("Configuration is empty, aborting")
  end

  for simulation in keys(benchmark_config)
    simulation_config = benchmark_config[simulation]

    for arch in keys(simulation_config)
      is_valid_arch = validate_arch(arch)
      if !is_valid_arch
        error("Configuration on '$arch' for '$simulation' is not valid, aborting")
      end

      simulation_entry = simulation_config[arch]
      if isempty(simulation_entry)
        error("Configuration for '$simulation' on '$arch' is defined but empty, aborting")
      end
    end
  end
end

function load_save_benchmark_data(benchmark_config)
  for simulation in keys(benchmark_config)

    mkpath(srcdir(simulation))
    simulation_config = benchmark_config[simulation]

    for arch in keys(simulation_config)
      simulation_entry = simulation_config[arch]

      new_simconfig = process_simulation_config(simulation_entry)

      open(get_config(simulation, arch), "w") do io
        TOML.print(io, new_simconfig)
      end
    end
  end
end

function validate_arch(arch::String)
  return arch in get_supported_arches()
end

function process_simulation_config(entry)
  received_params = dict_list(entry)

  params_list = Dict()
  add_meta_data!(params_list, received_params)
  add_task_data!(params_list, received_params)

  params_list
end

function add_meta_data!(params_list, received_params)
  meta_data = Dict("fields" => join(keys(first(received_params)), ","))
  push!(params_list, string(0) => meta_data)
end

function add_task_data!(params_list, received_params)
  for (idx, param) in enumerate(received_params)
    push!(params_list, string(idx) => param)
  end
end

function get_float_types(;just_real::Bool=true)
  just_real ? ["Float32", "Float64"] : ["Float32", "Float64", "ComplexF32", "ComplexF64"]
end