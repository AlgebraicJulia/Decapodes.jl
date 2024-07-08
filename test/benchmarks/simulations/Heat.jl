function setup_benchmark(config, ::Val{:Heat})
  Heat = @decapode begin
      U::Form0
      ∂ₜ(U) == 100 * Δ(U)
  end

  eval(gensim(Heat, code_target = config.code_target, stateeltype = config.float_type))
end

function create_mesh(config, ::Val{:Heat})
  s = triangulated_grid(100, 100, 1, 1, Point2{config.float_type})
  sd = EmbeddedDeltaDualComplex2D{Bool, config.float_type, Point2{config.float_type}}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(sd[:point]) do (x,_)
    return x
  end
  if config.code_target isa CUDATarget
    U = CuArray(U)
  end
  u₀ = ComponentArray(U=U)

  (sd, u₀)
end

function create_simulate(config, sd, simulate, ::Val{:Heat})
  simulate(sd, nothing, DiagonalHodge())
end

function run_simulation(config, fm, u0, ::Val{:Heat})
  tₑ = config.float_type(10)
  prob = ODEProblem(fm, u0, (0, tₑ), ())

  soln = solve(prob, Tsit5())
end
