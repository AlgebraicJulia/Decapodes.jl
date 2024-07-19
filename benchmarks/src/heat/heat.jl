using DrWatson
@quickactivate :benchmarks

struct HeatConfig
  float_type
  code_target
  res
end

function setup_config(task_config_data::Dict{String, Any})
  float_type = parse_float_type(task_config_data["float_type"])
  code_target = parse_code_target(task_config_data["code_target"])
  resolution = Float64(task_config_data["resolution"])
  @info "Float type: $(float_type), Code target: $(code_target), Resolution: $(resolution)"

  HeatConfig(float_type, code_target, resolution)
end

function setup_benchmark(config::HeatConfig)
  Heat = @decapode begin
      U::Form0
      ∂ₜ(U) == 100 * Δ(U)
  end

  eval(gensim(Heat, code_target = config.code_target, stateeltype = config.float_type))
end

function create_mesh(config::HeatConfig)
  s = triangulated_grid(100, 100, config.res, config.res, Point2{config.float_type})
  sd = EmbeddedDeltaDualComplex2D{Bool, config.float_type, Point2{config.float_type}}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(sd[:point]) do (x,_)
    return x
  end
  if config.code_target isa CUDATarget
    U = CuArray(U)
  end
  u₀ = ComponentArray(U=U)

  (sd, u₀, ())
end

function create_simulate(config::HeatConfig, sd, simulate)
  simulate(sd, nothing, DiagonalHodge())
end

function run_simulation(config::HeatConfig, fm, u0, cnst_param)
  tₑ = config.float_type(10)
  prob = ODEProblem(fm, u0, (0, tₑ), cnst_param)

  soln = solve(prob, Tsit5(), saveat=0.1)
end
