using DrWatson
@quickactivate :benchmarks

@info "Precompiling Julia"
using Pkg
Pkg.instantiate()
Pkg.precompile()
@info "Finished precompiling Julia"

using TOML

include(helpersdir("main_source.jl"))

if length(ARGS) == 0
    @info "Running all sims"
    generate_all_configs()
    run_all_physics()

elseif length(ARGS) == 3
    @info "Running single sim"
    generate_all_configs()

    const physics = ARGS[1]
    const arch = ARGS[2]
    const tag = ARGS[3]

    run_physics = SimNameData(physics, arch, tag)
    is_valid_config_instance(run_physics)
    run_single_physics(physics, [run_physics])
else
    error("Usage: ['physics' 'architecture' 'tag']")
end
