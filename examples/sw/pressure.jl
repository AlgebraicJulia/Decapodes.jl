using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
using Catlab.ACSetInterface
using CairoMakie
using Logging
using ComponentArrays

using GeometryBasics: Point3
Point3D = Point3{Float64}

PressureFlow = quote
  ## state variables
  V::Form1{X}
  P::Form0{X}

  ## derived quantities
  ΔV::Form1{X}
  ∇P::Form1{X}
  ΔP::Form0{X}
  ϕₚ::Form1{X}

  ## tanvars
  V̇::Form1{X}
  Ṗ::Form0{X}
  ∂ₜ(V) == V̇
  ∂ₜ(P) == Ṗ
  
  ∇P == d₀(P)
  ΔV == Δ₁(V)
  ΔP == Δ₀(P)

  V̇  == α(∇P) + μ(ΔV)
  ϕₚ == γ(-(L₀(V, P))) 
  Ṗ == β(Δ₀(P)) + ∘(dual_d₁,⋆₀⁻¹)(ϕₚ)
  ## Ṗ  == ϕₚ
end

pf = SummationDecapode(parse_decapode(PressureFlow))

flatten_form(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :k => x->2000x
    :μ => x->-0.0001x
    :α => x->0*x
    :β => x->2000*x
    :γ => x->1*x
    :i₀ => i0 
    :debug => (args...)->begin println(args[1], length.(args[2:end])) end
    _ => default_dec_generate(sd, my_symbol, hodge)
  end
  return (args...) ->  op(args...)
end

#include("coordinates.jl")

radius = 6371+90

primal_earth = loadmesh(Icosphere(3, radius))
nploc = argmax(x -> x[3], primal_earth[:point])

orient!(primal_earth)
primal_earth[:edge_orientation] = false
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())

physics = SummationDecapode(parse_decapode(PressureFlow))
gensim(expand_operators(physics), [:P, :V])
sim = eval(gensim(expand_operators(physics), [:P, :V]))

fₘ = sim(earth, generate)

begin
  vmag = 500
  velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4, 0))

  ## visualize the vector field
  ps = earth[:point]
  ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)
  arrows(
      ps, ns, fxaa=true, ## turn on anti-aliasing
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

u₀ = ComponentArrays(P = c, V=collect(v))
mesh(primal_earth, color=u₀.P, colormap=:plasma)
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
mass(soln, t, mesh, concentration=:P) = sum(⋆(0, mesh)*soln(t).concentration)

@show extrema(mass(soln, t, earth, :P) for t in 0:tₑ/150:tₑ)
end
mesh(primal_earth, color=soln(0).P, colormap=:jet)
mesh(primal_earth, color=soln(0) - soln(tₑ).P, colormap=:jet)
begin
## Plot the result
times = range(0.0, tₑ, length=150)
colors = [soln(t).P for t in times]

## Initial frame
fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = (-0.0001, 0.0001), colormap=:jet)
Colorbar(fig[1,2], ob)
framerate = 5

## Animation
record(fig, "weather.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = soln(t).P
end
end

AdvDiff = quote
    C::Form0{X}
    Ċ::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}
    starC::DualForm2{X}
    lvc::Form1{X}
    ## Fick's first law
    ϕ₁ ==  ∘(d₀,k,⋆₁)(C)
    ## Advective Flux
    ϕ₂ == -(L₀(V, C))
    ## Superposition Principle
    ϕ == plus(ϕ₁ , ϕ₂)
    ## Conservation of Mass
    Ċ == ∘(dual_d₁,⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

NavierStokes = quote
  V::Form1{X}
  V̇::Form1{X}
  G::Form1{X}
  T::Form0{X}
  ρ::Form0{X}
  ṗ::Form0{X}
  p::Form0{X}
  

  V̇ == neg₁(L₁′(V, V)) + 
        div₁(kᵥ(Δ₁(V) + third(d₀(δ₁(V)))), avg₀₁(ρ)) +
        d₀(half(i₁′(V, V))) +
        neg₁(div₁(d₀(p),avg₀₁(ρ))) +
        G
  ∂ₜ(V) == V̇
  ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))## + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
  ∂ₜ(p) == ṗ
end

parse_decapode(NavierStokes)
SummationDecapode(parse_decapode(NavierStokes))


