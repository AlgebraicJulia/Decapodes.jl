# Imports
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
using Logging
using JLD2
using Printf
using ComponentArrays
using GeometryBasics: Point3
Point3D = Point3{Float64}

# Define Model
Brusselator = @decapode begin
  (U, V)::Form0
  α::Constant
  F::Parameter

  U2V == (U*U) * V
  ∂ₜ(U)== 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  ∂ₜ(V) == (3.4 * U) - U2V + (α * Δ(U))
end
infer_types!(Brusselator)
resolve_overloads!(Brusselator)

# Define Domain
download("https://graphics.stanford.edu/courses/cs148-10-summer/as3/code/as3/teapot.obj", "teapot.obj")
s = EmbeddedDeltaSet2D("teapot.obj")
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())
fig,ax,ob = wireframe(sd)
save("teapot_subdivided.png", fig)

# Create initial data.
U = map(p -> abs(p[2]), point(s))
V = map(p -> abs(p[1]), point(s))
F₁ = map(sd[:point]) do (_,_,z)
  z ≥ 0.8 ? 5.0 : 0.0
end
F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₂ : F₁)
u₀ = ComponentArray(U=U, V=V)

# Run
function generate(sd, my_symbol; hodge=GeometricHodge()) end
sim = eval(gensim(Brusselator))
fₘ = sim(sd, generate)
tₑ = 1.15
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

@save "brusselator_teapot.jld2" soln

# Create side-by-side GIF
begin
frames = 800
# Initial frame
fig = CairoMakie.Figure(resolution = (1200, 1200))
p1 = CairoMakie.mesh(fig[1,1], s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
p2 = CairoMakie.mesh(fig[2,1], s, color=soln(0).V, colormap=:jet, colorrange=extrema(soln(0).V))
Colorbar(fig[1,2])
Colorbar(fig[2,2])

# Animation
record(fig, "brusselator_teapot.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    p1.plot.color = soln(t).U
    p2.plot.color = soln(t).V
end

end

