function setup_heat_benchmark(config)
  Heat = @decapode begin
      U::Form0
      ∂ₜ(U) == 100 * Δ(U)
  end

  eval(gensim(Heat, code_target = config.code_target, stateeltype = config.float_type))
end

function create_heat_mesh(config)
  s = triangulated_grid(100, 100, 1, 1, config.point_type)
  sd = EmbeddedDeltaDualComplex2D{Bool, config.float_type, config.point_type}(s)
  subdivide_duals!(sd, config.center_type)

  U = map(sd[:point]) do (x,_)
    return x
  end
  u₀ = ComponentArray(U=U)

  (sd, u₀)
end

function create_heat_simulate(config, sd, simulate)
  # hodge_type = @match center_type begin
  #   ::Circumcenter => DiagonalHodge()
  #   ::Barycenter => GeometricHodge()
  #   _ => error("Unsupported center type $center_type")
  # end
  simulate(sd, nothing, config.hodge_type)
end

function run_heat_simulation(config, fm, u0)
  tₑ = config.timeto_execute
  prob = ODEProblem(fm, u0, (0, tₑ), ())

  soln = solve(prob, Tsit5())
  soln.stats
end
