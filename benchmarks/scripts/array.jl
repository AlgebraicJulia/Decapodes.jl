using DrWatson
@quickactivate :benchmarks

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

const task_key = ARGS[1]
const sim_name = ARGS[2]
const arch = ARGS[3]

function extract_task_config(config_data)
  if !haskey(config_data, task_key)
    error("Warning: Task with key $(task_key) could not find config data")
  end

  task_config_data = all_config_data[task_key]
  @info string(task_config_data)

  task_config_data
end

function extract_stagenames(config_data)
  meta_config = all_config_data[string(0)]
  split(meta_config["stages"], ",")
end

@info "Running $sim_name on $arch, array id is $task_key"

# Extract data
all_config_data = TOML.parsefile(get_config(sim_name, arch))

task_config_data = extract_task_config(all_config_data)
solver_stages = get_solver_stages()

# Grab user's physics file
include(get_simfile(sim_name))

config = setup_config(task_config_data)

# Get intermediate variables to use in benchmarking
sim = setup_benchmark(config);
sd, u0, cnst_param = create_mesh(config);
fm = create_simulate(config, sd, sim);
result = run_simulation(config, fm, u0, cnst_param);

# Save solver statistics
stats_data = tostringdict(struct2dict(result.stats))
wsave(get_statsfile(task_key, sim_name, arch), stats_data)

# Setup and run benchmarking
simulation_suite = BenchmarkGroup()

simulation_suite[task_key][solver_stages[1]] = @benchmarkable setup_benchmark($config) gctrial=true
simulation_suite[task_key][solver_stages[2]] = @benchmarkable create_mesh($config) gcsample=true
simulation_suite[task_key][solver_stages[3]] = @benchmarkable create_simulate($config, $sd, $sim) gctrial=true
simulation_suite[task_key][solver_stages[4]] = @benchmarkable run_simulation($config, $fm, $u0, $cnst_param) gcsample=true

tune!(simulation_suite)
deca_sim_results = run(simulation_suite; verbose = true)
BenchmarkTools.save(get_benchfile(task_key, sim_name, arch), deca_sim_results)
