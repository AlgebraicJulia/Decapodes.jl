# Demonstrate how to encode boundary conditions in Decapodes using the "collage"
# technique.
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using GLMakie
using Logging
using JLD2
using Printf

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

BrusselatorDynamics = @decapode begin
  (U, V)::Form0
  U2V::Form0
  α::Constant
  F::Parameter
  U2V == (U .* U) .* V
  ∂ₜ(U) == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  ∂ₜ(V) == (3.4 * U) - U2V + (α * Δ(U))
end
# Visualize. You must have graphviz installed.
to_graphviz(BrusselatorDynamics)

# This is a "discrete" Decapode, with no morphisms.
BrusselatorBoundaries = @decapode begin
  L::Constant
end

# Specify the BC morphism between Decapodes.
BrusselatorBCMorphism = BCMorphism(ACSetTransformation(
  BrusselatorBoundaries, BrusselatorDynamics,
  Var = [1]))

# This is a "discrete" Decapode, with no morphisms.
BrusselatorInitialConditions = @decapode begin
  (U₀, V₀)::Form0
end

# Specify the IC morphism between Decapodes.
BrusselatorICMorphism = ICMorphism(ACSetTransformation(
  BrusselatorInitialConditions, BrusselatorDynamics,
  Var = [1,2]))

# Wrap these morphisms into a single collage.
BrusselatorCollage = Collage(
  BrusselatorBCMorphism,
  BrusselatorICMorphism)

# Create the BC loader, IC loader and simulation generator.
bc_loader, ic_loader, sim =
  simulation_helper(BrusselatorCollage)

# Specify a function that finds boundaries and values to assign, generic to any
# mesh.
function left_wall_constant_1(sd)
  min_x = minimum(x -> x[1], sd[:point])
  return (
    findall(x -> x[1] == min_x, sd[:point]),
    ones(count(x -> x[1] == min_x, sd[:point])))
end

# Specify a function that assigns initial values, generic to any mesh.
# This could be an explicit function, as below, or a function that
# interpolates data from e.g. a netcdf file.
function parabola_in_y(sd)
  map(sd[:point]) do (_,y)
    22 * (y *(1-y))^(3/2)
  end
end
function parabola_in_x(sd)
  map(sd[:point]) do (x,_)
    27 * (x *(1-x))^(3/2)
  end
end

# Specify some mesh.
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

# Demonstrate the boundary condition loader.
# This parameter data will be passed normally for demonstration.
F₁ = map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end
F₂ = zeros(nv(sd))

bc_generator = bc_loader(Dict(
    :L => left_wall_constant_1))

p = bc_generator(sd,
  (α = 0.001,
   F = t -> t ≥ 1.1 ? F₂ : F₁))

# Demonstrate the initial condition loader.
ic_generator = ic_loader(Dict(
  :U₀ => parabola_in_x,
  :V₀ => parabola_in_y))
u₀ = ic_generator(sd)

# Generate the simulation.
fₘ = eval(sim)(sd, default_dec_generate)

# Create problem and run sim for t ∈ [0,tₑ).
tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), p)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), p)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator_bounded.jld2" soln

# Visualize the final conditions.
GLMakie.mesh(s, color=findnode(soln(tₑ), :U), colormap=:jet)

begin # BEGIN Gif creation
frames = 100
# Initial frame
fig = GLMakie.Figure(resolution = (1200, 800))
p1 = GLMakie.mesh(fig[1,2], s, color=findnode(soln(0), :U), colormap=:jet, colorrange=extrema(findnode(soln(0), :U)))
p2 = GLMakie.mesh(fig[1,4], s, color=findnode(soln(0), :V), colormap=:jet, colorrange=extrema(findnode(soln(0), :V)))
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

# Animation
record(fig, "brusselator_bounded.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    p1.plot.color = findnode(soln(t), :U)
    p2.plot.color = findnode(soln(t), :V)
    lab1.text = @sprintf("%.2f",t)
end

end # END Gif creation
