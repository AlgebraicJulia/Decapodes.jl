#STOKES
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form
using ComponentArrays
using StaticArrays
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using ComponentArrays
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis

using Logging
using JLD2
using Printf

using GeometryBasics: Point2
Point2D = Point2{Float64}

StokesDynamics = @decapode begin
  (P)::Form0 ## Pressure.
  (v)::Form1 ## Velocity.
  (φ)::Constant
  (μ)::Constant

  ## Associate tangent variables with a state variable.
  ∂ₜ(v) == v̇
  ∂ₜ(P) == Ṗ
  
  v̇ == μ * Δ(v)-d₀(P) + φ
  Ṗ == ⋆₀⁻¹(dual_d₁(⋆₁(v)))
end 
# Visualize. You must have graphviz installed.
to_graphviz(StokesDynamics)
symsim = gensim(StokesDynamics)
s = triangulated_grid(1,1,1/100,1/100,Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())
f = evalsim(StokesDynamics)(sd,nothing)
g = (du,u) -> f(du,u,(φ=zeros(ne(s)),μ=1),0) 
ω = eval_constant_primal_form(sd,SVector{3}([1.,1.,1.]))
u = ComponentArray(v=ω,P=ones(nv(sd)))
du=similar(u)
g(du,u)
du

# This is a "discrete" Decapode, with no morphisms.
# TODO: Create an example with values that are not source terms.
BrusselatorBoundaries = @decapode begin
  Left::Constant
  Right::Constant
  Top::Constant
  Bottom::Constant
end

BrusselatorMorphism = @relation () begin
  rbl(C, Cl)
end

BrusselatorCollage = collate(
  BrusselatorDynamics,
  BrusselatorBoundaries,
  BrusselatorMorphism,
  Dict(
    :C => :U,
    :Cl => :L))
to_graphviz(BrusselatorCollage)

Brusselator = BrusselatorCollage

infer_types!(Brusselator)
resolve_overloads!(Brusselator)
to_graphviz(Brusselator)

s = triangulated_grid(1,1,1/10,1/10,Point3D)
#scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
#1/maximum(x->x[2], s[:point]),
#1.0])
#s[:point] = map(x -> scaling_mat*x, s[:point])
s[:edge_orientation] = false
orient!(s)

#Visualize the mesh.
wireframe(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())

#Select the vertices along the left "wall" of the domain.
diam(s) = sqrt(sum((s[:point][s[:∂v0][1]]- s[:point][s[:∂v1][1]]) .^2))
left_wall_idxs(s) = findall(x -> abs(x[1]) < diam(s), s[:point])
right_wall_idxs(s) = findall(x -> abs(x[1]-1) < diam(s), s[:point])
bottom_wall_idxs(s) = findall(x -> abs(x[2]) < diam(s), s[:point])
top_wall_idxs(s) = findall(x -> abs(x[2]-1) < diam(s), s[:point])
function generate(sd, my_symbol; hodge=GeometricHodge())
op = @match my_symbol begin
  #=
  :rbl => (u,u₀) -> begin
  u[left_wall_idxs(sd)] .= u₀
  u
  end
  :rbr => (u,u₀) -> begin
  u[right_wall_idxs(sd)] .= u₀
  u
  end
  :rbb => (u,u₀) -> begin
  u[bottom_wall_idxs(sd)] .= u₀
  u
  end
  :rbt => (u,u₀) -> begin
  u[top_wall_idx(sd)] .= u₀
  u
  end
  =#
  :.* => (x,y) -> x .* y
  _ => error("Unmatched operator $my_symbol")
end

return (args...) -> op(args...)
end



# Create initial data.
@assert all(map(sd[:point]) do (x,y)
  0.0 ≤ x ≤ 1.0 && 0.0 ≤ y ≤ 1.0
end)

U = map(sd[:point]) do (_,y)
  22 * (y *(2-y))^(3/2)
end

VV = map(sd[:point]) do (x,_)
  27 * (x *(2-x))^(3/2)
end

L = Float64[]#fill(1.0, length(left_wall_idxs(s)))

F₁ = map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end
mesh(s, color=F₁, colormap=:jet)

F₂ = zeros(nv(sd))

constants_and_parameters = (
  fourfour = 4.4,
  threefour = 3.4,
  α = 0.001,
  L = L,
  F = t -> t ≥ 1.1 ? F₂ : F₁)

# Generate the simulation.
simsym = gensim(StokesDynamics)
sim = eval(simsym)
fₘ = sim(sd, generate) #taking the dual mesh

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=VV)

# Visualize the initial conditions.
# If GLMakie throws errors, then update your graphics drivers,
# or use an alternative Makie backend like CairoMakie.
fig_ic = Figure()
p1 = mesh(fig_ic[1,2], s, color=u₀.U, colormap=:jet)
p2 = mesh(fig_ic[1,3], s, color=u₀.V, colormap=:jet)

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator_bounded.jld2" soln

# Visualize the final conditions.
mesh(s, color=soln(tₑ).U, colormap=:jet)

begin # BEGIN Gif creation
frames = 100
## Initial frame
fig = Figure(resolution = (1200, 800))
p1 = mesh(fig[1,2], s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
p2 = mesh(fig[1,4], s, color=soln(0).V, colormap=:jet, colorrange=extrema(soln(0).V))
ax1 = Axis(fig[1,2], width = 400, height = 400)
ax2 = Axis(fig[1,4], width = 400, height = 400)
hidedecorations!(ax1)
hidedecorations!(ax2)
hidespines!(ax1)
hidespines!(ax2)
Colorbar(fig[1,1])
Colorbar(fig[1,5])
Label(fig[1,2,Top()], "U")
Label(fig[1,4,Top()], "V")
lab1 = Label(fig[1,3], "")

## Animation
record(fig, "brusselator_bounded.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    p1.plot.color = soln(t).U
    p2.plot.color = soln(t).V
    lab1.text = @sprintf("%.2f",t)
end

end 
## END Gif creation

# Use boundaries on the left and right walls
BrusselatorBoundaries = @decapode begin
  L::Constant
  R::Constant
end

BrusselatorMorphism = @relation () begin
  rbl(C, Cl)
  rbr(C, Cr)
end

BrusselatorCollage = collate(
  BrusselatorDynamics,
  BrusselatorBoundaries,
  BrusselatorMorphism,
  Dict(
    :C => :U,
    :Cl => :L,
    :Cr => :R))
to_graphviz(BrusselatorCollage)

Brusselator = BrusselatorCollage

infer_types!(Brusselator)
resolve_overloads!(Brusselator)
to_graphviz(Brusselator)

gensim(expand_operators(Brusselator))
sim = eval(gensim(expand_operators(Brusselator)))
fₘ = sim(sd, generate)

L = fill(0.0, length(left_wall_idxs))
R = fill(0.0, length(right_wall_idxs))

constants_and_parameters = (
  fourfour = 4.4,
  threefour = 3.4,
  α = 0.001,
  L = L,
  R = R,
  F = t -> t ≥ 1.1 ? F₂ : F₁)

tₑ = 21.5e2

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator_double_bounded.jld2" soln

begin 
## BEGIN Gif creation
frames = 100
## Initial frame
fig = Figure(resolution = (1200, 800))
p1 = mesh(fig[1,2], s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
p2 = mesh(fig[1,4], s, color=soln(0).V, colormap=:jet, colorrange=extrema(soln(0).V))
ax1 = Axis(fig[1,2], width = 400, height = 400)
ax2 = Axis(fig[1,4], width = 400, height = 400)
hidedecorations!(ax1)
hidedecorations!(ax2)
hidespines!(ax1)
hidespines!(ax2)
Colorbar(fig[1,1])
Colorbar(fig[1,5])
Label(fig[1,2,Top()], "U")
Label(fig[1,4,Top()], "V")
lab1 = Label(fig[1,3], "")

## Animation
record(fig, "brusselator_dual_bounded.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    p1.plot.color = soln(t).U
    p2.plot.color = soln(t).V
    lab1.text = @sprintf("%.2f",t)
end

end 
## END Gif creation
