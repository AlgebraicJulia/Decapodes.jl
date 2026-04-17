using DrWatson
@quickactivate :benchmarks

using TOML

export mainsim_config_path, load_main_config, listof_main_physics,
  collect_simsfor_physics, has_config_args, config_args, config_arg

mainsim_config_path() = srcdir("main_config.toml")

load_main_config() = TOML.parsefile(mainsim_config_path())

listof_main_physics() = collect(keys(load_main_config()))

collect_simsfor_physics(physics) = collect_simsfor_physics(load_main_config(), physics)

function collect_simsfor_physics(main_config_info, physics)
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
  tag_info = physics_config_arch_info(physics_info, arch)
  for tag in keys(tag_info)
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

function has_config_args(config_info, snd::SimNameData)
  is_supported_arch(snd.arch) || return false

  haskey(config_info, snd.physics) || return false
  physics_info = config_info[snd.physics]

  haskey(physics_info, snd.arch) || return false
  arch_info = physics_info[snd.arch]

  haskey(arch_info, snd.tag) || return false
  tag_info = arch_info[snd.tag]

  return !isempty(tag_info)
end

function config_args(config_info, snd::SimNameData)
  if !(has_config_args(config_info, snd))
    error("Arguments for $(snd) were not found in the main configuration or provided architecture $(snd.arch) is invalid")
  end
  return config_info[snd.physics][snd.arch][snd.tag]
end

function config_arg(config_info, snd::SimNameData, arg)
  if !has_config_args(config_info, snd)
    return nothing
  end

  args = config_args(config_info, snd)
  if !haskey(args, arg)
    return nothing
  end

  return args[arg]
end
