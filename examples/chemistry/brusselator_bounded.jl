# TODO: Clean comments/ turn into a doc.
# Demonstrate how to encode boundary conditions in Decapodes using the "collage"
# technique.
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using ComponentArrays
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using MultiScaleArrays
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis

using Logging
using JLD2
using Printf

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

BrusselatorDynamics = @decapode begin
  ## Values living on vertices.
  (U, V)::Form0 ## State variables.
  U2V::Form0 ## Named intermediate variables.
  (U̇, V̇)::Form0 ## Tangent variables.
  (α)::Constant
  (F)::Parameter
  ## A named intermediate variable.
  U2V == (U .* U) .* V
  ## Specify how to compute the tangent variables.
  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(V))
  ## Associate tangent variables with a state variable.
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end
# Visualize. You must have graphviz installed.
to_graphviz(BrusselatorDynamics)

# This is a "discrete" Decapode, with no morphisms.
# TODO: Create an example with values that are not source terms.
BrusselatorBoundaries = @decapode begin
  L::Constant
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

# TODO: Create square domain of approximately 32x32 vertices.
s = triangulated_grid(30,10,2,2)
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0])
s[:point] = map(x -> scaling_mat*x, s[:point])
s[:edge_orientation] = false
orient!(s)
# Visualize the mesh.
wireframe(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())

# Select the vertices along the left "wall" of the domain.
min_x = minimum(x -> x[1], s[:point])
max_x = maximum(x -> x[1], s[:point])
left_wall_idxs = findall(x -> x[1] == min_x, s[:point])
right_wall_idxs = findall(x -> x[1] == max_x, s[:point])
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :rbl => (x,y) -> begin
      x[left_wall_idxs] .= y
      x
    end
    :rbr => (x,y) -> begin
      x[right_wall_idxs] .= y
      x
    end
    :.* => (x,y) -> x .* y
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

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

L = fill(1.0, length(left_wall_idxs))

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
gensim(expand_operators(Brusselator))
sim = eval(gensim(expand_operators(Brusselator)))
fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=V)

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
