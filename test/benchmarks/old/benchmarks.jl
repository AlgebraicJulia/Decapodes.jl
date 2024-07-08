module Benchmarks

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
using CUDA, CUDA.CUSPARSE

struct BenchConfig
  name::String

  float_type::DataType
  point_type
  center_type
  hodge_type

  mesh_size
  mesh_resolution

  code_target
  timeto_execute::AbstractFloat
end

# HOW TO ADD BENCHMARK
# 1. Setup should have decapode creation and run eval(gensim), return that func
# 2. Mesh creation should create the mesh and inital conditions, return both
# 3. Simulate creation should run the simulate from above with correct arguments]
# 4. Run should take the simulation code and initial conditions and run the solve,
# write the stats to a log
include("Heat.jl")
deca_sim_suite = BenchmarkGroup()

open("decapode_results_log.txt", "w") do f
  write(f, "Second to last number is count of accepted steps.\n")
end

for target in [CPUTarget(), CUDATarget()]
  for res in [1, 0.8, 0.5, 0.4, 0.2]
    config = BenchConfig("Heat", Float64, Point3, Circumcenter(), DiagonalHodge(), 100, res, target, 11.5)
    dispatch = Val(Symbol(config.name))

    sim = setup_benchmark(config, dispatch);
    sd, u0 = create_mesh(config, dispatch);
    fm = create_simulate(config, sd, sim, dispatch);
    result = run_simulation(config, fm, u0, dispatch)

    open("decapode_results_log.txt", "a") do f
      write(f, string(config)*"\n")
      write(f, string(result.stats)*"\n\n")
    end  
    
    deca_sim_suite[config.name][config.code_target]["Res: $(config.mesh_resolution)"]["Setup"] = @benchmarkable setup_benchmark($config, $dispatch) gctrial=true
    deca_sim_suite[config.name][config.code_target]["Res: $(config.mesh_resolution)"]["Mesh"] = @benchmarkable create_mesh($config, $dispatch) gcsample=true
    deca_sim_suite[config.name][config.code_target]["Res: $(config.mesh_resolution)"]["Simulate"] = @benchmarkable create_simulate($config, $sd, $sim, $dispatch) gctrial=true
    deca_sim_suite[config.name][config.code_target]["Res: $(config.mesh_resolution)"]["Solve"] = @benchmarkable run_simulation($config, $fm, $u0, $dispatch) gcsample=true
  end
end

@info "Tuning benchmarks"
tune!(deca_sim_suite)

@info "Saving parameters"
BenchmarkTools.save("decapode_bench_params.json", params(deca_sim_suite));

@info "Running benchmarks"
deca_sim_results = run(deca_sim_suite, verbose = true)
@info "Done running benchmarks"

@info "Saving results"
BenchmarkTools.save("decapode_results.json", deca_sim_results)
# BenchmarkTools.load("decapode_results.json")
end

# Gensim
# Simulate
# Mesh creation
# Solve

# By simulation
# By mesh size
# By processor
