# TODO: Do we have anything better solution than a runtime interface?

export setup_config, setup_benchmark, create_mesh, create_simulate, run_simulation

setup_config(args...) = begin
    error("Please add the following to your simulation file:'setup_config(task_config_data::Dict{String, Any})'")
end

setup_benchmark(args...) = begin
    error("Please add the following to your simulation file:'setup_benchmark(config::YourConfig)'")
end

create_mesh(args...) = begin
    error("Please add the following to your simulation file:'create_mesh(config::YourConfig)'")
end

create_simulate(args...) = begin
    error("Please add the following to your simulation file:'create_simulate(config::YourConfig, sd, simulate)'")
end

run_simulation(args...) = begin
    error("Please add the following to your simulation file: 'run_simulation(config::YourConfig, fm, u0, cnst_param)'")
end

