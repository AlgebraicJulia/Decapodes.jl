using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
using GLMakie
using Logging
# using CairoMakie 

using GeometryBasics: Point3
Point3D = Point3{Float64}

# Begin Navier Stokes

# Navier Stokes example
DiffusionExprBody = quote
  (T, Ṫ)::Form0{X}
  ϕ::DualForm1{X}

  # Fick's first law
  ϕ ==  ⋆₁(k(d₀(T)))
  # Diffusion equation
  Ṫ == ⋆₀⁻¹(dual_d₁(ϕ))
end
#DiffusionExprBody = quote
#  (T, Ṫ)::Form0{X}
#  ϕ::DualForm1{X}
#  k::Constant{X}
#
#  # Fick's first law
#  ϕ ==  ⋆₁(k*(d₀(T)))
#  # Diffusion equation
#  Ṫ == ⋆₀⁻¹(dual_d₁(ϕ))
#end

Diffusion = SummationDecapode(parse_decapode(DiffusionExprBody))
to_graphviz(Diffusion)

#AdvectionExprBody = quote
#  (T, Ṫ)::Form0{X}
#  V::Form1{X}
#  Ṫ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(T))))
#end
#
#Advection = SummationDecapode(parse_decapode(AdvectionExprBody))
#
#SuperpositionExprBody = quote
#  Ṫ₁::Form0{X}
#  Ṫ₂::Form0{X}
#  Ṫ::Form0{X}
#  T::Form0{X}
#  Ṫ == Ṫ₁ + Ṫ₂
#  ∂ₜ(T) == Ṫ
#end
#
#Superposition = SummationDecapode(parse_decapode(SuperpositionExprBody))

NavierStokesExprBody = quote
  (V, V̇, G)::Form1{X}
  (T, ρ, ṗ, p)::Form0{X}
  V̇ == neg₁(L₁′(V, V)) + 
        div₁(kᵥ(Δ₁(V) + third(d₀(δ₁(V)))), avg₀₁(ρ)) +
        d₀(half(i₁′(V, V))) +
        neg₁(div₁(d₀(p),avg₀₁(ρ))) +
        G
  ∂ₜ(V) == V̇
  ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
  ∂ₜ(p) == ṗ
end

NavierStokes = SummationDecapode(parse_decapode(NavierStokesExprBody))
NavierStokes[6, :type] = :Form0
NavierStokes[2, :type] = :Form1
to_graphviz(NavierStokes)
# TODO: Use infer_types!

EnergyExprBody = quote
  V::Form1{X}
  (ρ, p, T, Ṫ, Ṫₐ, Ṫ₁, bc₀)::Form0{X}

  ρ == div₀(p, R₀(T))
  Ṫₐ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(T))))
  bc₀ == ∂ₜₐ(Ṫₐ)
  Ṫ == Ṫₐ + Ṫ₁
  ∂ₜ(T) == Ṫ 
end

Energy = SummationDecapode(parse_decapode(EnergyExprBody))


# Needed until we resolve infered types
Energy[5, :type] = :Form0
to_graphviz(Energy)
# TODO: Use infer_types!

#BoundaryConditionsExprBody = quote
#  (V, V̇, bc₁)::Form1{X}
#  (Ṫ, ṗ, bc₀)::Form0{X}
#
#  # no-slip edges
#  bc₁ == ∂ᵥ(V̇)
#  # No change on left/right boundaries
#
#  # TODO: Change back to ∂ₜ once naming is fixed
#  bc₀ == ∂τ(Ṫ)
#  bc₀ == ∂ₚ(ṗ)
#end
#
#BoundaryConditions = SummationDecapode(parse_decapode(BoundaryConditionsExprBody))
#to_graphviz(BoundaryConditions)

compose_heat_xfer = @relation (V, ρ) begin
  flow(V, V̇, T, ρ, ṗ, p)
  energy(Ṫ, V, ρ, p, T, Ṫ₁)
  diffusion(T, Ṫ₁)
  #bcs(Ṫ, ṗ, V, V̇)
end
to_graphviz(compose_heat_xfer, junction_labels=:variable, box_labels=:name, prog="dot")

HeatXfer_comp = oapply(compose_heat_xfer,
                  [Open(NavierStokes, [:V, :V̇, :T, :ρ, :ṗ, :p]),
                   Open(Energy, [:Ṫ, :V, :ρ, :p, :T, :Ṫ₁]),
                   Open(Diffusion, [:T, :Ṫ])])
                   #Open(BoundaryConditions, [:Ṫ, :ṗ, :V, :V̇])])

HeatXfer = apex(HeatXfer_comp)
to_graphviz(HeatXfer)
to_graphviz(HeatXfer, directed=false)
# end

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :k => x->2000x
    :μ => x->-0.0001x
    # :μ => x->-2000x
    :α => x->0*x
    :β => x->2000*x
    :γ => x->1*x
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
    # :L₀ => (v,x)->dual_derivative(1, sd)*(i0(v, x))
    :L₀ => (v,x)->begin
      # dual_derivative(1, sd)*⋆(1, sd)*wedge_product(Tuple{1,0}, sd, v, x)
      ⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x)
    end
    :i₀ => i0 
    :Δ₀ => x -> begin # dδ
      δ(1, sd, d(0, sd)*x, hodge=hodge) end
    # :Δ₀ => x -> begin # d ⋆ d̃ ⋆⁻¹
    #   y = dual_derivative(1,sd)*⋆(1, sd, hodge=hodge)*d(0,sd)*x
    #   inv_hodge_star(0,sd, y; hodge=hodge)
    # end
    :Δ₁ => x -> begin # dδ + δd
      δ(2, sd, d(1, sd)*x, hodge=hodge) + d(0, sd)*δ(1, sd, x, hodge=hodge)
    end

    :δ₁ => x -> inv_hodge_star(0, sd, hodge=hodge) * dual_derivative(1,sd) * ⋆(1, sd, hodge=hodge) * x
    :i₁′ => (v,x) -> inv_hodge_star(0,sd, hodge=hodge) * wedge_product(Tuple{1,1}, sd, v, ⋆(1, sd, hodge=hodge) * x) #⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
    #:L₁′ = ... + d(0,sd)*i₁′(v,x) #i₀′(F1, d₁{X}(F1′)) + d₀{X}(i₁′(F1, F1′))
    :neg₁ => x -> -1.0 * x
    :neg₀ => x -> -1.0 * x
    :half => x -> 0.5 * x
    :third => x -> x / 3.0
    :div₀ => (v,x) -> v / x
    :div₁ => (v,x) -> v / x
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(s)
          append!(J, [s[e,:∂v0],s[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      sparse(I,J,V)*x
    end

    # :Δ₁ => x -> begin # d ⋆ d̃ ⋆⁻¹ + ⋆ d̃ ⋆ d
    #   y = dual_derivative(0,sd)*⋆(2, sd, hodge=hodge)*d(1,sd)*x
    #   inv_hodge_star(2,sd, y; hodge=hodge) 
    #   z = d(0, sd) * inverse_hode_star(2, sd, dual_derivative(1, sd)*⋆(1,sd, hodge=hodge)*x; hodge=hodge)
    #   return y + z
    # end
    :debug => (args...)->begin println(args[1], length.(args[2:end])) end
    x=> error("Unmatched operator $my_symbol")
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length.(args))"); op(args...);end
  return (args...) ->  op(args...)
end

# End Navier Stokes

flatten_form(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

#include("coordinates.jl")
#include("spherical_meshes.jl")

radius = 6371+90

#primal_earth = loadmesh(ThermoIcosphere())
primal_earth = loadmesh(Icosphere(3, radius))
nploc = argmax(x -> x[3], primal_earth[:point])

orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())

physics = SummationDecapode(parse_decapode(PressureFlow))
gensim(expand_operators(physics), [:P, :V])
sim = eval(gensim(expand_operators(physics), [:P, :V]))

fₘ = sim(earth, generate)

begin
  vmag = 500
  # velocity(p) = vmag*ϕhat(p)
  velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4, 0))
  # velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, -vmag/4, 0))

# visualize the vector field
  ps = earth[:point]
  ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)
  GLMakie.arrows(
      ps, ns, fxaa=true, # turn on anti-aliasing
      linecolor = :gray, arrowcolor = :gray,
      linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
      align = :center, axis=(type=Axis3,)
  )
end

begin
v = flatten_form(velocity, earth)
c_dist = MvNormal([radius/√(2), radius/√(2)], 20*[1, 1])
c = 100*[pdf(c_dist, [p[1], p[2]]) for p in earth[:point]]

theta_start = 45*pi/180
phi_start = 0*pi/180
x = radius*cos(phi_start)*sin(theta_start)
y = radius*sin(phi_start)*sin(theta_start)
z = radius*cos(theta_start)
c_dist₁ = MvNormal([x, y, z], 20*[1, 1, 1])
c_dist₂ = MvNormal([x, y, -z], 20*[1, 1, 1])

c_dist = MixtureModel([c_dist₁, c_dist₂], [0.6,0.4])

c = 100*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]


u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(collect(v))],Float64[], [:P, :V])
mesh(primal_earth, color=findnode(u₀, :P), colormap=:plasma)
tₑ = 30.0

@info("Precompiling Solver")
prob = ODEProblem(fₘ,u₀,(0,1e-4))
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())
@info("Done")
end

begin
mass(soln, t, mesh, concentration=:P) = sum(⋆(0, mesh)*findnode(soln(t), concentration))

@show extrema(mass(soln, t, earth, :P) for t in 0:tₑ/150:tₑ)
end
mesh(primal_earth, color=findnode(soln(0), :P), colormap=:jet)
mesh(primal_earth, color=findnode(soln(0) - soln(tₑ), :P), colormap=:jet)
begin
# Plot the result
times = range(0.0, tₑ, length=150)
colors = [findnode(soln(t), :P) for t in times]

# Initial frame
# fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = extrema(vcat(colors...)), colormap=:jet)
fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = (-0.0001, 0.0001), colormap=:jet)
Colorbar(fig[1,2], ob)
framerate = 5

# Animation
record(fig, "weather.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = findnode(soln(t), :P)
end
end

#AdvDiff = quote
#    C::Form0{X}
#    Ċ::Form0{X}
#    V::Form1{X}
#    ϕ::Form1{X}
#    ϕ₁::Form1{X}
#    ϕ₂::Form1{X}
#    starC::DualForm2{X}
#    lvc::Form1{X}
#    # Fick's first law
#    ϕ₁ ==  ∘(d₀,k,⋆₁)(C)
#    # Advective Flux
#    ϕ₂ == -(L₀(V, C))
#    # Superposition Principle
#    ϕ == plus(ϕ₁ , ϕ₂)
#    # Conservation of Mass
#    Ċ == ∘(dual_d₁,⋆₀⁻¹)(ϕ)
#    ∂ₜ(C) == Ċ
#end
#
#NavierStokes = quote
#  V::Form1{X}
#  V̇::Form1{X}
#  G::Form1{X}
#  T::Form0{X}
#  ρ::Form0{X}
#  ṗ::Form0{X}
#  p::Form0{X}
#  
#
#  V̇ == neg₁(L₁′(V, V)) + 
#        div₁(kᵥ(Δ₁(V) + third(d₀(δ₁(V)))), avg₀₁(ρ)) +
#        d₀(half(i₁′(V, V))) +
#        neg₁(div₁(d₀(p),avg₀₁(ρ))) +
#        G
#  ∂ₜ(V) == V̇
#  ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
#  ∂ₜ(p) == ṗ
#end
#
#parse_decapode(NavierStokes)
#SummationDecapode(parse_decapode(NavierStokes))

# Energy = @decapode Flow2DQuantities begin
#   (V)::Form1{X}
#   (ρ, p, T, Ṫ, Ṫₐ, Ṫ₁, bc₀)::Form0{X}

#   ρ == div₀(p, R₀(T))
#   Ṫₐ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(T))))
#   ∂ₜₐ(Ṫₐ) == bc₀
#   Ṫ == Ṫₐ + Ṫ₁
#   ∂ₜ{Form0{X}}(T) == Ṫ
# end

# BoundaryConditions = @decapode Flow2DQuantities begin 
#   (V, V̇, bc₁)::Form1{X}
#   (Ṫ, ṗ, bc₀)::Form0{X} 
#   # no-slip edges
#   ∂ᵥ(V̇) == bc₁
#   # No change on left/right boundaries
#   ∂ₜ(Ṫ) == bc₀
#   ∂ₚ(ṗ) == bc₀
# end
