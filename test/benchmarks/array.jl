using TOML
using BenchmarkTools
using ACSets
using Decapodes
using CombinatorialSpaces
using DiagrammaticEquations
using LinearAlgebra
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using GeometryBasics: Point2, Point3
# using CUDA, CUDA.CUSPARSE

arrayid = ARGS[1] # Job Array ID
sim_name = ARGS[2]
println(sim_name)
println(arrayid)
cd(sim_name)

println(pwd())

all_config_data = TOML.parsefile("$(sim_name)_cpu.toml")

task_key = string(arrayid)

if !haskey(all_config_data, task_key)
  error("Warning: Task with key $(task_key) could not find config data")
end

task_config_data = all_config_data[task_key]

println(string(task_config_data))

float_data = task_config_data["float_type"]
float_type = @match float_data begin
  "Float32" => Float32
  "Float64" => Float64
  "ComplexF32" => ComplexF32
  "ComplexF64" => ComplexF64
  _ => error("Float data $(float_data) is not valid, exiting early") 
end

code_target_data = task_config_data["code_target"]
code_target = @match code_target_data begin
  "CPUTarget" => CPUTarget()
  "CUDATarget" => CUDATarget()
  _ => error("Warning: Code target data $(code_target_data) is not valid, exiting early")
end

# TODO: This won't support Icospheres
resolution = parse(Float64, task_config_data["resolution"])

meta_config = all_config_data[string(0)]
solver_stages = split(meta_config["stages"], ",")

println("Float type: $(float_type), Code target: $(code_target), Resolution: $(resolution)")

struct BenchConfig
  float_type
  code_target
  res
end

config = BenchConfig(float_type, code_target, resolution)

include(joinpath(sim_name, "$sim_name.jl"))

simulation_suite = BenchmarkGroup()

sim = setup_benchmark(config);
sd, u0, cnst_param = create_mesh(config);
fm = create_simulate(config, sd, sim);
result = run_simulation(config, fm, u0, cnst_param);

remove_wrapper = r"^.*?\((.*?)\)$"
open(joinpath("results", "stats_$(task_key).txt"), "w") do file
  to_write = filter(!isspace, match(remove_wrapper, string(result.stats))[1])
  write(file, to_write)
end

simulation_suite[task_key][solver_stages[1]] = @benchmarkable setup_benchmark($config) gctrial=true
simulation_suite[task_key][solver_stages[2]] = @benchmarkable create_mesh($config) gcsample=true
simulation_suite[task_key][solver_stages[3]] = @benchmarkable create_simulate($config, $sd, $sim) gctrial=true
simulation_suite[task_key][solver_stages[4]] = @benchmarkable run_simulation($config, $fm, $u0, $cnst_param) gcsample=true

params_file_name = joinpath("params", "params_$(task_key).json")

if !isfile(params_file_name)
  println("Warning: Could not find previous parameters file, regenerating")
  tune!(simulation_suite)
  BenchmarkTools.save(params_file_name, params(simulation_suite))
end

deca_sim_results = run(simulation_suite, verbose = true)
BenchmarkTools.save(joinpath("results", "benchmarks_$(task_key).json"), deca_sim_results)
