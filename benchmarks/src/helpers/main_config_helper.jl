using DrWatson
@quickactivate :benchmarks

using TOML

export collect_mainconfig_sims, collect_mainconfig_simentries

function collect_mainconfig_sims(sim_name, main_config_info)
    entries = SimNameData[]
    haskey(main_config_info, sim_name) || return entries
    sim_config_info = main_config_info[sim_name]
  
    
    for arch in keys(sim_config_info)
        for tag in sim_config_info[arch]
           push!(entries, SimNameData(sim_name, arch, tag))
        end
    end

    return entries
end

function collect_mainconfig_simentries(sim_name, main_config_info)
  
    entries = SimNameData[]
    configured_sims = collect_mainconfig_sims(sim_name, main_config_info)
    for tagged_sim_namedata in configured_sims

        tagged_sim_config = config_path(tagged_sim_namedata)
        if !isfile(tagged_sim_config) 
            @info "Config file for $sim_name on $arch named $tag not found, skipping"
            continue
        end
    
        tagged_sim_data = TOML.parsefile(tagged_sim_config)
        num_entries = get_autoconfig_size(tagged_sim_data)
        for task_id in 1:num_entries
            
            sim_name = tagged_sim_namedata.sim_name
            arch = tagged_sim_namedata.arch
            tag = tagged_sim_namedata.tag
            task_key = string(task_id)

            push!(entries, SimNameData(sim_name, arch, tag, task_key))
        end
    end
  
    entries
end
  