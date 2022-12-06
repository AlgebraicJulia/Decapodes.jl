using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
#using GLMakie
using Logging
using CairoMakie 

using GeometryBasics: Point3
Point3D = Point3{Float64}

PressureFlow = quote
  # state variables
  V::Form1{X}
  P::Form0{X}
  C::Form0{X}

  # parameters
  k::Constant{ℝ}
  μ::Constant{ℝ}
  α::Constant{ℝ}
  βc::Constant{ℝ}
  γc::Constant{ℝ}
  βₚ::Constant{ℝ}
  γₚ::Constant{ℝ}


  # derived quantities
  ΔV::Form1{X}
  ∇P::Form1{X}
  ∇C::Form1{X}
  ΔP::Form0{X}
  ΔC::Form0{X}
  ϕₚ::Form1{X}
  ϕc::Form1{X}
  # tanvars
  V̇::Form1{X}
  Ṗ::Form0{X}
  Ċ::Form0{X}
  ∂ₜ(V) == V̇
  ∂ₜ(P) == Ṗ
  ∂ₜ(C) == Ċ
  
  ∇P == d₀(P)
  ∇C == d₀(C)
  ΔV == Δ₁(V)
  ΔP == Δ₀(P)
  ΔC == Δ₀(C)

  V̇  == α*∇P + μ*ΔV
  ϕₚ == γₚ*(-(L₀(V, P)))
  ϕc == γc*(-(L₀(V, C))) 
  Ṗ == βₚ*Δ₀(P) + ∘(dual_d₁,⋆₀⁻¹)(ϕₚ)
  Ċ == βc*Δ₀(C) + ∘(dual_d₁,⋆₀⁻¹)(ϕc)
  
  
  # Ṗ  == ϕₚ
end

pf = SummationDecapode(parse_decapode(PressureFlow))

flatten_form(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :⋆₀ => x->⋆(0,sd,hodge=hodge)*x
    :⋆₁ => x->⋆(1, sd, hodge=hodge)*x
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=hodge)
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=hodge)*x
    :d₀ => x->d(0,sd)*x
    :d₁ => x->d(1,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
    :plus => (+)
    :(-) => x-> -x
    :L₀ => (v,x)->begin
      ⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x)
    end
    :i₀ => i0 
    :Δ₀ => x -> begin # dδ
      δ(1, sd, d(0, sd)*x, hodge=hodge) end
    :Δ₁ => x -> begin # dδ + δd
      δ(2, sd, d(1, sd)*x, hodge=hodge) + d(0, sd)*δ(1, sd, x, hodge=hodge)
    end

    :debug => (args...)->begin println(args[1], length.(args[2:end])) end
    x=> error("Unmatched operator $my_symbol")
  end
  return (args...) ->  op(args...)
end


radius = 6371+90

primal_earth = loadmesh(Icosphere(4,radius))
nploc = argmax(x -> x[3], primal_earth[:point])

orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3}(primal_earth)
subdivide_duals!(earth, Circumcenter())

physics = SummationDecapode(parse_decapode(PressureFlow))
gensim(expand_operators(physics))
sim = eval(gensim(expand_operators(physics)))

fₘ = sim(earth, generate)

vmag = 5
velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4, 0))


begin
  v = flatten_form(velocity, earth)

  theta_start = 45*pi/180
  phi_start = 45*pi/180
  x = radius*cos(phi_start)*sin(theta_start)
  y = radius*sin(phi_start)*sin(theta_start)
  z = radius*cos(theta_start)

  #c_dist₁ = MvNormal([x, y, z], 200*[1, 1, 1])
  #c_dist₂ = MvNormal([x, y, -z], 200*[1, 1, 1])

  #c_dist = MixtureModel([c_dist₁, c_dist₂], [0.6,0.4])
  #c = 1e20*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]

  c_dist = MvNormal([0, 0, 0], 20*[radius, radius, radius/40])
  c = 1e20*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]

  p_dist = MvNormal([0, 0, 0], 20*[radius, radius, radius/40])
  pfield = 1e20*[pdf(p_dist, [p[1], p[2], p[3]]) for p in earth[:point]]
  #pfield = 100000*[p[3]/radius for p in earth[:point]]
  #maximum(c)
  extrema(pfield)
  fig, ax, ob = mesh(primal_earth, color=pfield)
  Colorbar(fig[1,2], ob)
  CairoMakie.save("p_field.png", fig)


  u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(collect(v)), VectorForm(pfield)],Float64[], [:C, :V, :P])
  #fig, ax, ob = mesh(primal_earth, color=findnode(u₀, :P), colormap=:plasma)
  #Colorbar(fig[1,2], ob)
end
extrema(pfield)

d(0, earth, pfield)

tₑ = 10.0

for α_prime in 1000
  var_name = "alpha_"
  println(α_prime)
  my_consts = (k=1,μ=-1,α=α_prime,βc=1e5,γc=-1,βₚ=1e5,γₚ=-1)
  #check stability
  prob = ODEProblem(fₘ,u₀,(0,1e-4),my_consts)
  soln = solve(prob, Tsit5())
  soln.retcode != :Unstable || error("Solver was not stable")
  #solve the actual ODE
  prob = ODEProblem(fₘ,u₀,(0,tₑ),my_consts)
  soln = solve(prob, Tsit5())
  soln.retcode != :Unstable || error("Solver was not stable")

  numframes = 1000;

  # Plot the result
  times = range(0.0, tₑ, length=numframes)
  colors = [findnode(soln(t), :C) for t in times]

  for i in eachindex(colors)
    if any(isnan.(colors[i]))
      println("found NaN")
      continue
    end
  end


  # Initial frame
  #fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = extrema(vcat(colors...)), colormap=:jet)
  fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = (5e4, 2e5), colormap=:jet)
  Colorbar(fig[1,2], ob)
  #framerate = 5

  # Animation
  record(fig, "weather_" *var_name*string(α_prime)* "_11_29_b_1e5_gc_m1_gp_m1_equator_0p.gif", range(0.0, tₑ; length=numframes); framerate = 120) do t
      ob.color = findnode(soln(t), :C)
    end

end 



