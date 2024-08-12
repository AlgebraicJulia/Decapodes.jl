using DrWatson
@quickactivate :benchmarks
using TOML

function process_benchmark_config(simconfig_path)
    benchmark_config = TOML.parsefile(simconfig_path)
    validate_config(benchmark_config)
    load_save_benchmark_data(benchmark_config)
end

# TODO: Update this to deal with different tagged sims
# Keep validation seperate to prevent creation of intermediate files/folders if error
function validate_config(benchmark_config)
  if isempty(benchmark_config)
    error("Configuration is empty")
  end

  for simulation in keys(benchmark_config)
    simulation_config = benchmark_config[simulation]

    for arch in keys(simulation_config)
      is_valid_arch = is_supported_arch(arch)
      if !is_valid_arch
        error("Configuration on '$arch' for '$simulation' is not valid")
      end

      simulation_entry = simulation_config[arch]
      if isempty(simulation_entry)
        error("Configuration for '$simulation' on '$arch' is defined but empty")
      end
    end
  end
end

function load_save_benchmark_data(benchmark_config)
  for physics in keys(benchmark_config)

    physics_sim_entry = benchmark_config[physics]

    for arch in keys(physics_sim_entry)
      arch_sim_entry = physics_sim_entry[arch]

      for tag in keys(arch_sim_entry)
        tagged_sim_entry = arch_sim_entry[tag]
        new_simconfig = process_simulation_config(tagged_sim_entry)

        sim_namedata = SimNameData(physics, arch, tag)
        open(simconfig_path(sim_namedata), "w") do io
          TOML.print(io, new_simconfig)
        end
      end

    end
  end
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
