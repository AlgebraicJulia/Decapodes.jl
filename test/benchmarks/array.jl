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
println(arrayid)
println(pwd())

all_config_data = TOML.parsefile("benchmarks_config_heat.toml")

task_key = string(arrayid)

if !haskey(all_config_data, task_key)
  print("Warning: Task with key $(task_key) could not find config data, exiting early")
  exit()
end

task_config_data = all_config_data[task_key]

println(string(task_config_data))

float_data = task_config_data["float_type"]
float_type = @match float_data begin
  "Float32" => Float32
  "Float64" => Float64
  _ => println("Warning: Float data $(float_data) is not valid, exiting early")
end

code_target_data = task_config_data["code_target"]
code_target = @match code_target_data begin
  "CPUTarget" => CPUTarget()
  "CUDATarget" => CUDATarget()
  _ => println("Warning: Code target data $(code_target_data) is not valid, exiting early")
end

println("Float type: $(float_type), Code target type: $(code_target)")

struct BenchConfig
  name::String
  float_type::DataType
  code_target
end

config = BenchConfig("Heat", float_type, code_target)

include(joinpath("simulations", "Heat.jl"))

simulation_suite = BenchmarkGroup()

dispatch = Val(Symbol(config.name))

sim = setup_benchmark(config, dispatch);
sd, u0 = create_mesh(config, dispatch);
fm = create_simulate(config, sd, sim, dispatch);
# result = run_simulation(config, fm, u0, dispatch);

simulation_suite[task_key]["Setup"] = @benchmarkable setup_benchmark($config, $dispatch) gctrial=true
simulation_suite[task_key]["Mesh"] = @benchmarkable create_mesh($config, $dispatch) gcsample=true
simulation_suite[task_key]["Simulate"] = @benchmarkable create_simulate($config, $sd, $sim, $dispatch) gctrial=true
simulation_suite[task_key]["Solve"] = @benchmarkable run_simulation($config, $fm, $u0, $dispatch) gcsample=true

params_file_name = joinpath("params", "benchmark_heat_params_$(task_key).json")

if !isfile(params_file_name)
  println("Warning: Could not find previous parameters file, regenerating")
  tune!(simulation_suite)
  BenchmarkTools.save(params_file_name, params(simulation_suite))
end

deca_sim_results = run(simulation_suite, verbose = true)
BenchmarkTools.save(joinpath("results", "benchmark_heat_results_$(task_key).json"), deca_sim_results)
