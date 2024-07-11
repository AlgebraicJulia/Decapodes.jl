function setup_benchmark(config)
  Heat = @decapode begin
      U::Form0
      ∂ₜ(U) == 100 * Δ(U)
  end

  eval(gensim(Heat, code_target = config.code_target, stateeltype = config.float_type))
end

function create_mesh(config)
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

function create_simulate(config, sd, simulate)
  simulate(sd, nothing, DiagonalHodge())
end

function run_simulation(config, fm, u0, cnst_param)
  tₑ = config.float_type(10)
  prob = ODEProblem(fm, u0, (0, tₑ), cnst_param)

  soln = solve(prob, Tsit5())
end
