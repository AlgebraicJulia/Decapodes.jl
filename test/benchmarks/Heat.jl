function setup_benchmark(config, ::Val{:Heat})
  Heat = @decapode begin
      U::Form0
      ∂ₜ(U) == 100 * Δ(U)
  end

  eval(gensim(Heat, code_target = config.code_target, stateeltype = config.float_type))
end

function create_mesh(config, ::Val{:Heat})
  s = triangulated_grid(config.mesh_size, config.mesh_size, config.mesh_resolution, config.mesh_resolution, config.point_type{config.float_type})
  sd = EmbeddedDeltaDualComplex2D{Bool, config.float_type, config.point_type{config.float_type}}(s)
  subdivide_duals!(sd, config.center_type)

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
  # hodge_type = @match center_type begin
  #   ::Circumcenter => DiagonalHodge()
  #   ::Barycenter => GeometricHodge()
  #   _ => error("Unsupported center type $center_type")
  # end
  simulate(sd, nothing, config.hodge_type)
end

function run_simulation(config, fm, u0, ::Val{:Heat})
  tₑ = config.float_type(config.timeto_execute)
  prob = ODEProblem(fm, u0, (0, tₑ), ())

  soln = solve(prob, Tsit5())
  
  open("decapode_results_log.txt", "a") do f
    write(f, string(config)*"\n")
    write(f, string(soln.stats)*"\n\n")
  end

end
