export solver_stages, supported_arches, is_supported_arch, meta_config_id, physics_config_filename

const solver_stages_list = ["Setup", "Mesh", "Simulate", "Solve"]
solver_stages() = return solver_stages_list

const supported_architectures_list = ["cpu", "cuda"]
supported_arches() = return supported_architectures_list
is_supported_arch(arch) = return arch in supported_arches()

const meta_key = string(0)
meta_config_id() = return meta_key

const physconfig_filename = "config.toml"
physics_config_filename() = return physconfig_filename
