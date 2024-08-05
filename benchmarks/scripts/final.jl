using DrWatson
@quickactivate :benchmarks

include(helpersdir("data_aggr_helper.jl"))

using TOML

const slurm_id = ARGS[1]
const sim_name = ARGS[2]

aggregate_data(slurm_id, sim_name)

run(`julia --threads=auto $(postprocessdir("default_out.jl")) $slurm_id $sim_name`)