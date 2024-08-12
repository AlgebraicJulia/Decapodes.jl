using DrWatson
@quickactivate :benchmarks

using TOML

function collect_entriesfor_physics(main_config_info, physics)
  entries = SimNameData[]
  add_entriesfor_physics!(entries, main_config_info, physics)
end

function add_entriesfor_physics!(entries, main_config_info, physics)
  physics_info = main_config_physics_info(main_config_info, physics)
  for arch in keys(physics_info)
    add_entriesfor_physics_arch!(entries, physics_info, physics, arch)
  end
  return entries
end

function add_entriesfor_physics_arch!(entries, physics_info, physics, arch)
  physics_arch_taglist = physics_config_arch_info(physics_info, arch)
  for tag in physics_arch_taglist
    add_entryfor_tagged!(entries, physics, arch, tag)
  end
  return entries
end

function add_entryfor_tagged!(entries, physics, arches, tag)
  push!(entries, SimNameData(physics, arches, tag))
  return entries
end

function main_config_physics_info(main_config_info, physics)
  haskey(main_config_info, physics) || error("Physics $physics does not exist in the provided main configuration")
  return main_config_info[physics]
end

function physics_config_arch_info(physics_info, arch)
  haskey(physics_info, arch) || error("Architecture $arch does not exist in the provided physics entry")
  physics_arch_taglist = physics_info[arch]
  return physics_arch_taglist
end
