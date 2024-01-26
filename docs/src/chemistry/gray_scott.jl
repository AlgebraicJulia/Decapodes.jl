using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using GLMakie
using Logging
using JLD2
using Printf
using ComponentArrays

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

# We use the model equations as stated here and use the initial conditions for
# f, k, rᵤ, rᵥ as listed for experiment 4.
# https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/
GrayScott = @decapode begin
  (U, V)::Form0
  (UV2)::Form0
  (U̇, V̇)::Form0
  (f, k, rᵤ, rᵥ)::Constant
  UV2 == (U .* (V .* V))
  U̇ == rᵤ * Δ(U) - UV2 + f * (1 .- U)
  V̇ == rᵥ * Δ(V) + UV2 - (f + k) .* V
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

# Visualize. You must have graphviz installed.
to_graphviz(GrayScott)

# We resolve types of intermediate variables using sets of rules.
infer_types!(GrayScott)
to_graphviz(GrayScott)

# Resolve overloads. i.e. ~dispatch
resolve_overloads!(GrayScott)
to_graphviz(GrayScott)

s = loadmesh(Rectangle_30x10())
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0])
s[:point] = map(x -> scaling_mat*x, s[:point])
s[:edge_orientation] = false
orient!(s)
# Visualize the mesh.
GLMakie.wireframe(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
function generate(sd, my_symbol; hodge=GeometricHodge()) end

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

constants_and_parameters = (
  f = 0.024,
  k = 0.055,
  rᵤ = 0.01,
  rᵥ = 0.005)

# Generate the simulation.
gensim(expand_operators(GrayScott))
sim = eval(gensim(expand_operators(GrayScott)))
fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U,V=V)

# Visualize the initial conditions.
# If GLMakie throws errors, then update your graphics drivers,
# or use an alternative Makie backend like CairoMakie.
fig_ic = GLMakie.Figure()
p1 = GLMakie.mesh(fig_ic[1,2], s, color=u₀.U, colormap=:jet)
p2 = GLMakie.mesh(fig_ic[1,3], s, color=u₀.V, colormap=:jet)
display(fig_ic)

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "gray_scott.jld2" soln

# Visualize the final conditions.
GLMakie.mesh(s, color=soln(tₑ).U, colormap=:jet)

begin # BEGIN Gif creation
frames = 100
# Initial frame
fig = GLMakie.Figure(resolution = (1200, 800))
p1 = GLMakie.mesh(fig[1,2], s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
p2 = GLMakie.mesh(fig[1,4], s, color=soln(0).V, colormap=:jet, colorrange=extrema(soln(0).V))
ax1 = Axis(fig[1,2], width = 400, height = 400)
ax2 = Axis(fig[1,4], width = 400, height = 400)
hidedecorations!(ax1)
hidedecorations!(ax2)
hidespines!(ax1)
hidespines!(ax2)
Colorbar(fig[1,1], colormap=:jet, colorrange=extrema(soln(0).U))
Colorbar(fig[1,5], colormap=:jet, colorrange=extrema(soln(0).V))
Label(fig[1,2,Top()], "U")
Label(fig[1,4,Top()], "V")
lab1 = Label(fig[1,3], "")

# Animation
record(fig, "gray_scott.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    p1.plot.color = soln(t).U
    p2.plot.color = soln(t).V
    lab1.text = @sprintf("%.2f",t)
end

end # END Gif creation

# Run on the sphere.
# You can use lower resolution meshes, such as Icosphere(3).
s = loadmesh(Icosphere(5))
orient!(s)
# Visualize the mesh.
GLMakie.wireframe(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
function generate(sd, my_symbol; hodge=GeometricHodge()) end

# Create initial data.
U = map(sd[:point]) do (_,y,_)
  abs(y)
end

V = map(sd[:point]) do (x,_,_)
  abs(x)
end

constants_and_parameters = (
  f = 0.024,
  k = 0.055,
  rᵤ = 0.01,
  rᵥ = 0.005)

# Generate the simulation.
fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U,V=V)

# Visualize the initial conditions.
# If GLMakie throws errors, then update your graphics drivers,
# or use an alternative Makie backend like CairoMakie.
fig_ic = GLMakie.Figure()
p1 = GLMakie.mesh(fig_ic[1,2], s, color=u₀.U, colormap=:jet)
p2 = GLMakie.mesh(fig_ic[1,3], s, color=u₀.V, colormap=:jet)
display(fig_ic)

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "gray_scott_sphere.jld2" soln

# Visualize the final conditions.
GLMakie.mesh(s, color=soln(tₑ).U, colormap=:jet)

begin # BEGIN Gif creation
frames = 800
# Initial frame
fig = GLMakie.Figure(resolution = (1200, 1200))
p1 = GLMakie.mesh(fig[1,1], s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
p2 = GLMakie.mesh(fig[2,1], s, color=soln(0).V, colormap=:jet, colorrange=extrema(soln(0).V))
Colorbar(fig[1,2], colormap=:jet, colorrange=extrema(soln(0).U))
Colorbar(fig[2,2], colormap=:jet, colorrange=extrema(soln(0).V))

# Animation
record(fig, "gray_scott_sphere.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    p1.plot.color = soln(t).U
    p2.plot.color = soln(t).V
end

end # END Gif creation
