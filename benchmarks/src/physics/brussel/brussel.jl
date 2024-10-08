using DrWatson
@quickactivate :benchmarks

using ACSets
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

using Decapodes.Canon.Chemistry

function pass_simulation_instance()
  return SimulationInstance(setup_config, setup_benchmark, create_mesh, create_simulate, run_simulation)
end

struct BrusselConfig
  float_type
  code_target
  res
end

function setup_config(task_config_data::Dict{String, Any})
  float_type = parse_float_type(task_config_data["float_type"])
  code_target = parse_code_target(task_config_data["code_target"])
  resolution = Float64(task_config_data["resolution"])
  @info "Float type: $(float_type), Code target: $(code_target), Resolution: $(resolution)"

  BrusselConfig(float_type, code_target, resolution)
end

function setup_benchmark(config::BrusselConfig)
  eval(gensim(brusselator, code_target = config.code_target, stateeltype = config.float_type))
end

function create_mesh(config::BrusselConfig)
  s = triangulated_grid(1, 1, config.res, config.res, Point2{config.float_type})
  sd = EmbeddedDeltaDualComplex2D{Bool, config.float_type, Point2{config.float_type}}(s)
  subdivide_duals!(sd, Circumcenter())

  U = config.float_type.(map(sd[:point]) do (_,y)
    22 * (y *(1-y))^(3/2)
  end)

  V = config.float_type.(map(sd[:point]) do (x,_)
    27 * (x *(1-x))^(3/2)
  end)

  F₁ = config.float_type.(map(sd[:point]) do (x,y)
    (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
  end)

  F₂ = zeros(config.float_type, nv(sd))

  if config.code_target isa CUDATarget
    U = CuArray(U)
    V = CuArray(V)
    F₁ = CuArray(F₁)
    F₂ = CuArray(F₂)
  end

  u₀ = ComponentArray(U=U, V=V)

  constants_and_parameters = (
    α = config.float_type.(0.001),
    F = t -> t ≥ config.float_type.(1.1) ? F₁ : F₂)

  (sd, u₀, constants_and_parameters)
end

function create_simulate(config::BrusselConfig, sd, simulate)
  simulate(sd, nothing, DiagonalHodge())
end

function run_simulation(config::BrusselConfig, fm, u0, cnst_param)
  tₑ = config.float_type(11.5)
  prob = ODEProblem(fm, u0, (0, tₑ), cnst_param)

  soln = solve(prob, Tsit5(), saveat=0.1)
end
