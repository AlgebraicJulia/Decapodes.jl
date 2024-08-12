using DrWatson
@quickactivate :benchmarks

using BenchmarkTools
using TOML

const task_key = ARGS[1]
const physics = ARGS[2]
const arch = ARGS[3]
const tag = ARGS[4]

sim_namedata = SimNameData(physics, arch, tag, task_key)

function extract_task_config(config_data)
  if !haskey(config_data, task_key)
    error("Warning: Task with key $(task_key) could not find config data")
  end

  task_config_data = all_config_data[task_key]
  @info string(task_config_data)

  task_config_data
end

@info "Running $physics on $arch, tagged as $tag, array id is $task_key"

# Extract data
all_config_data = TOML.parsefile(simconfig_path(sim_namedata))
task_config_data = extract_task_config(all_config_data)

# Grab user's physics file
include(physicspath(physics))

sim_instance = pass_simulation_instance()

config = sim_instance.setup_config(task_config_data)

# Get intermediate variables to use in benchmarking
sim = sim_instance.setup_benchmark(config);
sd, u0, cnst_param = sim_instance.create_mesh(config);
fm = sim_instance.create_simulate(config, sd, sim);
result = sim_instance.run_simulation(config, fm, u0, cnst_param);

# Save solver statistics
stats_data = tostringdict(struct2dict(result.stats))
wsave(statsfile_path(sim_namedata), stats_data)

# Setup and run benchmarking
simulation_suite = BenchmarkGroup()

stages = solver_stages()
simulation_suite[task_key][stages[1]] = @benchmarkable sim_instance.setup_benchmark($config) gctrial=true
simulation_suite[task_key][stages[2]] = @benchmarkable sim_instance.create_mesh($config) gcsample=true
simulation_suite[task_key][stages[3]] = @benchmarkable sim_instance.create_simulate($config, $sd, $sim) gctrial=true
simulation_suite[task_key][stages[4]] = @benchmarkable sim_instance.run_simulation($config, $fm, $u0, $cnst_param) gcsample=true

tune!(simulation_suite)
deca_sim_results = run(simulation_suite; verbose = true)
BenchmarkTools.save(benchfile_path(sim_namedata), deca_sim_results)
