function process_benchmark_config(configname)
    benchmark_config = TOML.parsefile(configname)
    validate_config(benchmark_config)
    load_benchmark_data(benchmark_config)
end

function load_benchmark_data(benchmark_config)
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

function validate_arch(arch::String)
  return arch == "cpu" || arch == "cuda"
end

function process_simulation_config(entry)
  received_params = dict_list(entry)

  params_list = Dict()

  meta_data = Dict("fields" => join(keys(first(received_params)), ","))
  push!(params_list, string(0) => meta_data)

  for (idx, param) in enumerate(received_params)
    push!(params_list, string(idx) => param)
  end

  params_list
end

function get_float_types(;just_real=true)
  just_real ? ["Float32", "Float64"] : ["Float32", "Float64", "ComplexF32", "ComplexF64"]
end
