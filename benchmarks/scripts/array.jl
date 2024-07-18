using DrWatson
@quickactivate "benchmarks"

using ACSets
using BenchmarkTools
using CUDA
using CUDA.CUSPARSE
using CombinatorialSpaces
using ComponentArrays
using Decapodes
using DiagrammaticEquations
using GeometryBasics: Point2, Point3
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using TOML

arrayid = ARGS[1] # Job Array ID
sim_name = ARGS[2]
code_target_name = ARGS[3]

println(sim_name)
println(arrayid)
println(code_target_name)

simdir(args...) = srcdir(sim_name, args...)
paramsdir(args...) = simdir("params", args...)
resultsdir(args...) = datadir("sims", sim_name, args...)

all_config_data = TOML.parsefile(simdir("$(sim_name)_$(code_target_name).toml"))

task_key = string(arrayid)

if !haskey(all_config_data, task_key)
  error("Warning: Task with key $(task_key) could not find config data")
end

task_config_data = all_config_data[task_key]
println(string(task_config_data))

meta_config = all_config_data[string(0)]
solver_stages = split(meta_config["stages"], ",")

include(simdir("$sim_name.jl"))

config = setup_config(task_config_data)

simulation_suite = BenchmarkGroup()

sim = setup_benchmark(config);
sd, u0, cnst_param = create_mesh(config);
fm = create_simulate(config, sd, sim);
result = run_simulation(config, fm, u0, cnst_param);

stats_data = tostringdict(struct2dict(result.stats))
wsave(resultsdir("stats_$(task_key)_$(code_target_name).jld2"), stats_data)

simulation_suite[task_key][solver_stages[1]] = @benchmarkable setup_benchmark($config) gctrial=true
simulation_suite[task_key][solver_stages[2]] = @benchmarkable create_mesh($config) gcsample=true
simulation_suite[task_key][solver_stages[3]] = @benchmarkable create_simulate($config, $sd, $sim) gctrial=true
simulation_suite[task_key][solver_stages[4]] = @benchmarkable run_simulation($config, $fm, $u0, $cnst_param) gcsample=true

# params_file_name = paramsdir("params_$(task_key)_$(code_target_name).json")

# if !isfile(params_file_name)
#   println("Warning: Could not find previous parameters file, regenerating")
#   tune!(simulation_suite)
#   BenchmarkTools.save(params_file_name, params(simulation_suite))
# end

# simulation_params = only(BenchmarkTools.load(params_file_name))

tune!(simulation_suite)
deca_sim_results = run(simulation_suite, verbose = true)
BenchmarkTools.save(resultsdir("benchmarks_$(task_key)_$(code_target_name).json"), deca_sim_results)
