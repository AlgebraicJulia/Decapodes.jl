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
using ComponentArrays
using GeometryBasics: Point3

Point3D = Point3{Float64}

flatten(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :k => x->2000x
    _ => default_dec_generate(sd, my_symbol, hodge)
  end

  return (args...) ->  op(args...)
end


DiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}

    ## Fick's first law
    ϕ ==  ∘(d₀, k)(C)
    ## Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end


diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
gensim(expand_operators(ddp), [:C])
f = eval(gensim(expand_operators(ddp), [:C]))

#include("coordinates.jl")

const RADIUS = 6371+90
#primal_earth, npi, spi = makeSphere(0, 180, 5, 0, 360, 5, RADIUS);
#nploc = primal_earth[npi, :point]
#primal_earth = loadmesh(ThermoIcosphere())
primal_earth = loadmesh(Icosphere(4, RADIUS))
nploc = argmax(x -> x[3], primal_earth[:point])
primal_earth[:edge_orientation] = false
orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())


fₘ = f(earth, generate)
c_dist = MvNormal(nploc[[1,2]], 100[1, 1])
c = [pdf(c_dist, [p[1], p[2]]./√RADIUS) for p in earth[:point]]

u₀ = ComponentArray(C=c)
tₑ = 10
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

using CairoMakie 
import CairoMakie: wireframe, mesh, Figure, Axis, arrows

#mesh(primal_earth, color=soln(0).C, colormap=:plasma)
#mesh(primal_earth, color=soln(tₑ).C, colormap=:plasma)
#mesh(primal_earth, color=soln(tₑ)-soln(0).C, colormap=:plasma)

begin
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

advdiff = parse_decapode(AdvDiff)
advdiffdp = SummationDecapode(advdiff)
gensim(expand_operators(advdiffdp), [:C, :V])
sim = eval(gensim(expand_operators(advdiffdp), [:C, :V]))

fₘ = sim(earth, generate)
end

begin
  vmag = 500
  velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4))

  ## visualize the vector field
  ps = earth[:point]
  ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)

  GLMakie.arrows(
      ps, ns, fxaa=true, ## turn on anti-aliasing
      linecolor = :gray, arrowcolor = :gray,
      linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
      align = :center, axis=(type=Axis3,)
  )
end

begin
v = flatten(velocity, earth)
c_dist = MvNormal([RADIUS/√(2), RADIUS/√(2)], 20*[1, 1])
c = 100*[pdf(c_dist, [p[1], p[2]]) for p in earth[:point]]

theta_start = 45*pi/180
phi_start = 0*pi/180
x = RADIUS*cos(phi_start)*sin(theta_start)
y = RADIUS*sin(phi_start)*sin(theta_start)
z = RADIUS*cos(theta_start)
c_dist₁ = MvNormal([x, y, z], 20*[1, 1, 1])
c_dist₂ = MvNormal([x, y, -z], 20*[1, 1, 1])

c_dist = MixtureModel([c_dist₁, c_dist₂], [0.6,0.4])

c = 100*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]

u₀ = ComponentArray(C=c,V=collect(v))
tₑ = 30.0

prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())
end

begin
mass(soln, t, mesh, concentration=:C) = sum(⋆(0, mesh)*soln(t).concentration)
@show extrema(mass(soln, t, earth, :C) for t in 0:tₑ/150:tₑ)
end

begin

## Plot the result
times = range(0.0, tₑ, length=150)
colors = [soln(t).C for t in times]

## Initial frame
fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = (-0.0001, 0.0001), colormap=:jet)
Colorbar(fig[1,2], ob)
framerate = 5

## Animation
record(fig, "diff_adv.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = soln(t).C
end
end
