using DrWatson
@quickactivate :benchmarks

include(helpersdir("data_aggr_helper.jl"))

using TOML

const slurm_id = ARGS[1]
const sim_name = ARGS[2]
const arch = ARGS[3]
const tag = ARGS[4]

aggregate_data(slurm_id, sim_name)

run(`julia --threads=auto $(postprocessdir("default_out.jl")) $slurm_id $sim_name`)