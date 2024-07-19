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

all_config_data = TOML.parsefile(get_config(sim_name, arch))

task_config_data = extract_task_config(all_config_data)
solver_stages = extract_stagenames(all_config_data)

include(get_simfile(sim_name))

config = setup_config(task_config_data)

sim = setup_benchmark(config);
sd, u0, cnst_param = create_mesh(config);
fm = create_simulate(config, sd, sim);
result = run_simulation(config, fm, u0, cnst_param);

stats_data = tostringdict(struct2dict(result.stats))
wsave(get_statsfile(task_key, sim_name, arch), stats_data)

simulation_suite = BenchmarkGroup()

simulation_suite[task_key][solver_stages[1]] = @benchmarkable setup_benchmark($config) gctrial=true
simulation_suite[task_key][solver_stages[2]] = @benchmarkable create_mesh($config) gcsample=true
simulation_suite[task_key][solver_stages[3]] = @benchmarkable create_simulate($config, $sd, $sim) gctrial=true
simulation_suite[task_key][solver_stages[4]] = @benchmarkable run_simulation($config, $fm, $u0, $cnst_param) gcsample=true

# params_file_name = paramsdir("params_$(task_key)_$(arch).json")

# if !isfile(params_file_name)
#   @info "Warning: Could not find previous parameters file, regenerating"
#   tune!(simulation_suite)
#   BenchmarkTools.save(params_file_name, params(simulation_suite))
# end

# simulation_params = only(BenchmarkTools.load(params_file_name))

tune!(simulation_suite)
deca_sim_results = run(simulation_suite, verbose = true)
BenchmarkTools.save(get_benchfile(task_key, sim_name, arch), deca_sim_results)
