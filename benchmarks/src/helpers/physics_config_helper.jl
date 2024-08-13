using DrWatson
@quickactivate :benchmarks
using TOML

function generate_configsfor_physics(physics)
    simulations = collect_simsfor_physics(physics)

    for sim in simulations
      if !has_config_args(load_physicsconfig(sim), sim)
        error("Simulation $(sim) does not have a valid configuration entry")
      end
    end

    for sim in simulations
      load_save_benchmark_data(sim)
    end
end

function load_save_benchmark_data(sim)
  task_simconfig = process_simulation_config(access_config_args(load_physicsconfig(sim), sim))
  open(simconfig_path(sim), "w") do io
    TOML.print(io, task_simconfig)
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
