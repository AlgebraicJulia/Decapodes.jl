using DrWatson
@quickactivate :benchmarks

include(helpersdir("data_aggr_helper.jl"))

using TOML

const slurm_id = ARGS[1]
const physics = ARGS[2]

aggregate_data(slurm_id, physics)

run(`julia --threads=auto $(postprocessdir("docs_pp.jl")) $slurm_id $physics`)
