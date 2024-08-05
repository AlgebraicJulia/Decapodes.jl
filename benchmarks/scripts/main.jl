using DrWatson
@quickactivate :benchmarks

@info "Precompiling Julia"
using Pkg
Pkg.instantiate()
Pkg.precompile()
@info "Finished precompiling Julia"

using TOML

function run_all_physics()
    physics_list = readdir(physicsdir())
    main_config_info = TOML.parsefile(mainconfig_path())
    validate_all_physics(physics_list, main_config_info)

    for physics in physics_list
        physics_configs = collect_mainconfig_sims(physics, main_config_info)
        run_single_physics(physics, physics_configs)
    end
end

function validate_all_physics(all_physics, main_config_info)
    for physics in all_physics
        physics_configs = collect_mainconfig_sims(physics, main_config_info)
        for config in physics_configs
            validate_config_instance(config)
        end
    end
end

function validate_config_instance(sim_namedata::SimNameData)
    sim_name = sim_namedata.sim_name

    is_supported_arch(sim_namedata.arch)  || error("Architecture $(arch) is not in list $(supported_arches())")
    simfile = get_simfile(sim_name)
    if !isfile(simfile)
        error("Simulation file at $(simfile) was not found")
    end

    config = config_path(sim_namedata)
    if !isfile(config)
        error("Config file at $(config) was not found")
    end
    return true
end

function run_single_physics(physics, physics_configs)
    
    if isempty(physics_configs)
        return
    end

    rm(resultsdir(physics), recursive=true, force=true)
    mkpath(resultsdir(physics))

    dependency_ids = []

    for config in physics_configs
        run_single_config!(dependency_ids, config)
    end

    if isempty(dependency_ids)
        return
    end

    run(`sbatch --dependency=afterok:$(join(dependency_ids, ",")) $(scriptsdir("final.sh")) $physics`)
end

function run_single_config!(dependency_ids, sim_namedata::SimNameData)
    config_data = TOML.parsefile(config_path(sim_namedata))

    concurrent_jobs = 6

    sim_name = sim_namedata.sim_name
    arch = sim_namedata.arch
    tag = sim_namedata.tag

    count = get_autoconfig_size(config_data)

    @info "Running $count $arch tasks"

    if arch == "cpu"
        jobid = readchomp(`sbatch --array=1-$(count)%$(concurrent_jobs) --parsable --partition=hpg-milan $(scriptsdir("array.sh")) $sim_name $arch $tag`)
    elseif arch == "cuda"
        jobid = readchomp(`sbatch --array=1-$(count)%$(concurrent_jobs) --parsable --partition=gpu --gres=gpu:a100:1 $(scriptsdir("array.sh")) $sim_name $arch $tag`)
    end

    @info "Job ID for $(tagfor_run(sim_namedata)) is: $jobid"
    push!(dependency_ids, jobid)
end

if length(ARGS) == 0
    @info "Running all sims"
    run_all_physics()
elseif length(ARGS) == 3
    @info "Running single sim"
    const sim_name = ARGS[1]
    const arch = ARGS[2]
    const tag = ARGS[3]

    run_physics = SimNameData(sim_name, arch, tag)
    validate_config_instance(run_physics)
    run_single_physics(sim_name, [run_physics])
else
    error("Usage: ['sim_name' 'architecture' 'tag']")
end

# @info "Regenerating configuration files"
# run(`julia --threads=auto config_generate.jl`)


