using DrWatson
@quickactivate "benchmarks"

using TOML
using MLStyle

const sim_name = ARGS[1]
const architecture = ARGS[2]
const stages = ["Setup", "Mesh", "Simulate", "Solve"]
const tracker = "clean.txt"

function get_code_target()
    @match architecture begin
        "cpu" => "CPUTarget"
        "cuda" => "CUDATarget"
        _ => error("Second argument should be either 'cpu' or 'cuda'")
    end
end

function get_float_types(;just_real=true)
    just_real ? ["Float32", "Float64"] : ["Float32", "Float64", "ComplexF32", "ComplexF64"]
end

function update_tracker(file_name)
    rm(file_name, force=true)
end

simdir(args...) = srcdir(sim_name, args...)

dictparams(_) = error("You have not yet implemented a parameter dictionary for this simulation")

function dictparams(::Val{:heat})
    resolution = [5, 2, 1]
    code_target = get_code_target()
    float_type = get_float_types()
    @strdict resolution code_target float_type
end

function generate_config(filename)
    mkpath(simdir())
    params = dict_list(dictparams(Val(Symbol(sim_name))))

    params_list = Dict()

    meta_data = Dict("fields" => join(keys(first(params)), ","), "stages" => join(stages, ","))
    push!(params_list, string(0) => meta_data)

    for (idx, param) in enumerate(params)
        push!(params_list, string(idx) => param)
    end

    open(joinpath(simdir(), filename), "w") do io
        TOML.print(io, params_list)
    end

    update_tracker(simdir(tracker))
end

generate_config("$(sim_name)_$(architecture).toml")


