using DrWatson
@quickactivate :benchmarks

using TOML
using MLStyle

const default_config_file = srcdir("main_config.toml")
# TODO: Simulation is too hard-coded to allow other stages right now
const stages = ["Setup", "Mesh", "Simulate", "Solve"]
const tracker = "clean.txt"

function get_float_types(;just_real=true)
    just_real ? ["Float32", "Float64"] : ["Float32", "Float64", "ComplexF32", "ComplexF64"]
end

function update_tracker(file_name)
    rm(file_name, force=true)
end

function process_benchmark_config(configname)
    benchmark_config = TOML.parsefile(configname)

    for simulation in keys(benchmark_config)
        mkpath(srcdir(simulation))
        simulation_config = benchmark_config[simulation]

        for arch in keys(simulation_config)
            simulation_entry = simulation_config[arch]

            if isempty(simulation_entry)
                error("Configuration for '$simulation' on '$arch' is defined but empty")
            end
            new_simconfig = process_simulation_config(simulation_entry)

            new_configname = get_configname(simulation, arch)
            open(srcdir(simulation, new_configname), "w") do io
                TOML.print(io, new_simconfig)
            end

            # update_tracker(simdir(tracker))
        end
    end
    
end

function process_simulation_config(entry)
    received_params = dict_list(entry)

    params_list = Dict()

    # TODO: Remove this maybe, really just needed for final table
    meta_data = Dict("fields" => join(keys(first(received_params)), ","), "stages" => join(stages, ","))
    push!(params_list, string(0) => meta_data)

    for (idx, param) in enumerate(received_params)
        push!(params_list, string(idx) => param)
    end

    params_list
end

if length(ARGS) == 0
    config_file = default_config_file
elseif length(ARGS) == 1
    config_file = ARGS[1]
else
    error("Usage: [config file name]")
end

process_benchmark_config(config_file)


