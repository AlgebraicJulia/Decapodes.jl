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

struct Config
  float_type::DataType
  point_type::DataType
  center_type

  hodge_type

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

heat_config = Config(Float64, Point3{Float64}, Circumcenter(), DiagonalHodge(), CPUTarget(), 11.5)

deca_sim_suite = BenchmarkGroup()

sim = setup_heat_benchmark(heat_config)
sd, u0 = create_heat_mesh(heat_config)
fm = create_heat_simulate(heat_config, sd, sim)
res = run_heat_simulation(heat_config, fm, u0)

deca_sim_suite["Heat"]["Setup"] = @benchmarkable setup_heat_benchmark(heat_config)
deca_sim_suite["Heat"]["Mesh"] = @benchmarkable create_heat_mesh(heat_config)
deca_sim_suite["Heat"]["Simulate"] = @benchmarkable create_heat_simulate(heat_config, sd, sim)
deca_sim_suite["Heat"]["Solve"] = @benchmarkable run_heat_simulation(heat_config, fm, u0)

@info "Simulation Benchmarks "

deca_sim_results = run(deca_sim_suite, verbose = true, seconds = 1)

for sim in sort(collect(keys(deca_sim_results)))
    test = median(deca_sim_results[sim])

    println("Simulation: $sim")
    for k in sort(collect(keys(test)))
        t = test[k].time / 1e6
        m = test[k].memory / 1e6
        println("Portion: $k, [$t ms, $m MB]")
    end
    println("----------------------------------------------------------------")
end

end

# Gensim
# Simulate
# Mesh creation
# Solve

# By simulation
# By mesh size
# By processor
