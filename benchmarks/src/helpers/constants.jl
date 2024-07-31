export get_solver_stages, get_supported_arches, meta_config_id

const solver_stages = ["Setup", "Mesh", "Simulate", "Solve"]
get_solver_stages() = return solver_stages

const supported_architectures = ["cpu", "cuda"]
get_supported_arches() = return supported_architectures

const meta_key = string(0)
meta_config_id() = return meta_key
