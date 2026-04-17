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

function pass_simulation_instance()
  return SimulationInstance(setup_config, setup_benchmark, create_mesh, create_simulate, run_simulation)
end

struct CHConfig
  float_type
  code_target
  res
end

function setup_config(task_config_data::Dict{String, Any})
  float_type = parse_float_type(task_config_data["float_type"])
  code_target = parse_code_target(task_config_data["code_target"])
  resolution = float_type(task_config_data["resolution"])
  @info "Float type: $(float_type), Code target: $(code_target), Resolution: $(resolution)"

  CHConfig(float_type, code_target, resolution)
end

function setup_benchmark(config::CHConfig)
  CahnHilliard = @decapode begin
    C::Form0
    (D, γ)::Constant
    ∂ₜ(C) == D * Δ(C.^3 - C - γ * Δ(C))
  end

  eval(gensim(CahnHilliard, code_target = config.code_target, stateeltype = config.float_type))
end

function create_mesh(config::CHConfig)
  s = triangulated_grid(100, 100, config.res, config.res, Point2{config.float_type})
  sd = EmbeddedDeltaDualComplex2D{Bool, config.float_type, Point2{config.float_type}}(s)
  subdivide_duals!(sd, Circumcenter())

  C = rand(config.float_type, nv(sd))

  if config.code_target isa CUDATarget
    C = CuArray(C)
  end
  u₀ = ComponentArray(C=C)

  constants = (D = config.float_type(0.5), γ = config.float_type(0.5))

  (sd, u₀, constants)
end

function create_simulate(config::CHConfig, sd, simulate)
  simulate(sd, nothing, DiagonalHodge())
end

function run_simulation(config::CHConfig, fm, u0, cnst_param)
  tₑ = config.float_type(200)
  prob = ODEProblem(fm, u0, (0, tₑ), cnst_param)

  soln = solve(prob, Tsit5(), saveat=0.1)
end
