using Catlab
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using JLD2
using Printf
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis
using ComponentArrays
using Symbolics
using LinearSolve

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

# Values living on vertices.
# State variables.
# Named intermediate variables.
Brusselator = @decapode begin
  (U, V)::Form0
  U2V::Form0
  (U̇, V̇)::Form0

  (α)::Constant
  F::Parameter

  U2V == (U .* U) .* V

  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(V))
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

# Visualize. You must have graphviz installed.
to_graphviz(Brusselator)

# We resolve types of intermediate variables using sets of rules.

infer_types!(Brusselator)
# Visualize. Note that variables now all have types.
to_graphviz(Brusselator)

# Resolve overloads. i.e. ~dispatch
resolve_overloads!(Brusselator)
# Visualize. Note that functions are renamed.
to_graphviz(Brusselator)

s = triangulated_grid(1,1,0.008,0.008,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());

# Create initial data.
@assert all(map(sd[:point]) do (x,y)
  0.0 ≤ x ≤ 1.0 && 0.0 ≤ y ≤ 1.0
end)

U = map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end

V = map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end

F₁ = map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end

F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₁ : F₂)

# Generate the simulation.
sim = evalsim(Brusselator)
fₘ = sim(sd, nothing, DiagonalHodge())

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=V)

tₑ = 11.5

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator.jld2" soln

# Visualize the final conditions.
fig = Figure();
ax = Axis(fig[1,1])
mesh!(ax, s, color=soln(tₑ).U, colormap=:jet)

function save_dynamics(save_file_name)
  time = Observable(0.0)
  u = @lift(soln($time).U)
  f = Figure()
  ax = CairoMakie.Axis(f[1,1], title = @lift("Brusselator U Concentration at Time $($time)"))
  gmsh = mesh!(ax, s, color=u, colormap=:jet,
               colorrange=extrema(soln(tₑ).U))
  Colorbar(f[1,2], gmsh)
  timestamps = range(0, tₑ, step=1e-1)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("brusselator_explicit.gif")

# This begins the demonstration of using implicit solvers
sim = evalsim(Brusselator, preallocate=false)
fₘ = sim(sd, nothing, DiagonalHodge())

u₀ = ComponentArray(U=U, V=V)

du₀ = copy(u₀)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> fₘ(du, u, constants_and_parameters, 0.0), du₀, u₀)

f = ODEFunction(fₘ; jac_prototype = float.(jac_sparsity))
tₑ = 11.5
@info("Solving")
prob = ODEProblem(f, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, FBDF(linsolve = KLUFactorization()), progress=true, progress_steps=1, saveat=0.1);
@info("Done")

save_dynamics("brusselator_implicit.gif")
# This ends the demonstration of using implicit solvers

# Run on the sphere.
# You can use lower resolution meshes, such as Icosphere(3).
s = loadmesh(Icosphere(7));
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

# Create initial data.
U = map(sd[:point]) do (_,y,_)
  abs(y)
end

V = map(sd[:point]) do (x,_,_)
  abs(x)
end

F₁ = map(sd[:point]) do (_,_,z)
  z ≥ 0.8 ? 5.0 : 0.0
end

F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₁ : F₂)

# Generate the simulation.
fₘ = sim(sd, nothing, DiagonalHodge())

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=V)

# Visualize the initial conditions.
fig = Figure()
ax = LScene(fig[1,1], scenekw=(lights=[],))
mesh!(ax, s, color=u₀.U, colormap=:jet)
display(fig)

tₑ = 11.5

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator_sphere.jld2" soln

# Visualize the final conditions.
fig = Figure()
ax = LScene(fig[1,1], scenekw=(lights=[],))
mesh!(ax, s, color=soln(tₑ).U, colormap=:jet)
display(fig)

function save_dynamics(save_file_name)
  time = Observable(0.0)
  u = @lift(soln($time).U)
  f = Figure()
  ax = LScene(f[1,1], scenekw=(lights=[],))
  gmsh = mesh!(ax, s, color=u, colormap=:jet,
               colorrange=extrema(soln(tₑ).U))
  Colorbar(f[1,2], gmsh)
  timestamps = range(0, tₑ, step=1e-1)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics("brusselator_sphere_U.gif")
