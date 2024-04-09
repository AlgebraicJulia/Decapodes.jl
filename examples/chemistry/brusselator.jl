using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using CairoMakie
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

Brusselator = @decapode begin
  # Values living on vertices.
  (U, V)::Form0{X} # State variables.
  (U2V)::Form0{X} # Named intermediate variables.
  (U̇, V̇)::Form0{X} # Tangent variables.
  # Scalars.
  (α)::Constant{X}
  F::Parameter{X}
  # A named intermediate variable.
  U2V == (U .* U) .* V
  # Specify how to compute the tangent variables.
  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(U))
  # Associate tangent variables with a state variable.
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

# TODO: Create square domain of approximately 32x32 vertices.
# s = loadmesh(Rectangle_30x10())
# Visualize the mesh.
s = triangulated_grid(200, 200, 1, 1, Point3D);
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0]);
s[:point] = map(x -> scaling_mat*x, s[:point]);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

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
  F = t -> t ≥ 1.1 ? F₂ : F₁)

gensim(Brusselator)
sim = evalsim(Brusselator)
fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=V)

tₑ = 25

begin
  @info("Precompiling Solver")
  prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
  soln = solve(prob, Tsit5())
  soln.retcode != :Unstable || error("Solver was not stable")
  @info("Solving")
  prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
  soln = solve(prob, Tsit5())
  @info("Done")
end

begin
    frames = 300
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    msh = CairoMakie.mesh!(ax, s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
    Colorbar(fig[1,2], msh)
    CairoMakie.record(fig, "Brusselator_CPU.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
        msh.color = soln(t).U
    end
end

# Run on the sphere.
# You can use lower resolution meshes, such as Icosphere(3).
s = loadmesh(Icosphere(5));

# Visualize the mesh.
wireframe(s)

sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter())

# Create initial data.
U = map(sd[:point]) do (_,y,_)
  abs(y)
end

V = map(sd[:point]) do (x,_,_)
  abs(x)
end

# TODO: Try making this sparse.
F₁ = map(sd[:point]) do (_,_,z)
  z ≥ 0.8 ? 5.0 : 0.0
end
mesh(s, color=F₁, colormap=:jet)

# TODO: Try making this sparse.
F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₁ : F₂
  )

# Generate the simulation.
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

# Visualize the final conditions.
mesh(s, color=soln(tₑ).U, colormap=:jet)

# BEGIN Gif creation
begin 
frames = 800
fig = Figure(resolution = (1200, 1200))
p1 = mesh(fig[1,1], s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
p2 = mesh(fig[2,1], s, color=soln(0).V, colormap=:jet, colorrange=extrema(soln(0).V))
Colorbar(fig[1,2])
Colorbar(fig[2,2])

record(fig, "brusselator_sphere.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    p1.plot.color = soln(t).U
    p2.plot.color = soln(t).V
end

end 
# END Gif creation
