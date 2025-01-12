# Brusselator

This Brusselator example is adapted from [MethodOfLines.jl](https://docs.sciml.ai/MethodOfLines/stable/tutorials/brusselator/#brusselator)'s page on the same topic. The [Brusselator](https://en.wikipedia.org/wiki/Brusselator) is a autocatalytic chemical reaction that takes place between two reactants `U` and `V`.

```@setup INFO
include(joinpath(Base.@__DIR__, ".." , "..", "docinfo.jl"))
info = DocInfo.Info()
```
## Dependencies
```@example DEC
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis

using Catlab
using CombinatorialSpaces
using ComponentArrays
using DiagrammaticEquations
using Decapodes
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}
nothing # hide
```

## The Model

We establish the model for the Brusselator, with the two reactants `U` and `V`, modeled as residing on the vertices of the mesh. The equations encode a reaction that occurs independently at each point coupled with a diffusion term as well as a source term `F` in the case of `U`. Here `α` denotes the rate of diffusion for both reactants.
```@example DEC
BrusselatorDynamics = @decapode begin
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
nothing # hide
```

## Boundary Conditions

We now establish the Dirichlet boundary conditions for our model. Here we intend to set some portion of the `U` variable to be a fixed value on some portion of the mesh. At this point the boundary conditions are only set symbolically and their actual implementation can change. Note that these values are set at the beginning of execution, as shown by the computation graph.
```@example DEC
BrusselatorBoundaries = @decapode begin
  B::Constant
end

BrusselatorMorphism = @relation () begin
  rlb(C, Cb)
end

Brusselator = collate(
  BrusselatorDynamics,
  BrusselatorBoundaries,
  BrusselatorMorphism,
  Dict(
    :C => :U,
    :Cb => :B))

to_graphviz(Brusselator)
```

## The Mesh

We load our triangulated mesh with horizontal and vertical resolution being `h=0.01`. `Point3D` is being used for the primal mesh `s` for ease of visualization while `Point2D` is used for the dual mesh `sd` for better memory usage. Since this conversion only drops the final z-coordinate, no information is lost.
```@example DEC
h = 0.01
s = triangulated_grid(1,1,h,h,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());

fig = Figure() # hide
ax = CairoMakie.Axis(fig[1,1], aspect=1) # hide
wf = wireframe!(ax, s; linewidth=1) # hide
save("Brusselator_rect.png", fig) # hide
nothing # hide
```

!["BrusselatorRect"](Brusselator_rect.png)

## Initial data
We assign the initial values of `U` and `V` according to a continuous function. Since they both live on the vertices of our mesh, we can simply iterate over all point coordinates, extract the coordinate values (for either x or y) and compute the desired value. `F` has some logic attached to it encoding that it will "activate" only once the simulation has reached time `t = 1.1`.

Here we also decide to set our boundary conditions to be `1.0` along the left and right sides of our mesh.
```@example DEC
U = map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end

V = map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end

fig = Figure() # hide
ax = CairoMakie.Axis(fig[1,1], aspect=1, title = "Initial value of U") # hide
msh = CairoMakie.mesh!(ax, s, color=U, colormap=:jet, colorrange=extrema(U)) # hide
Colorbar(fig[1,2], msh) # hide
save("initial_U.png", fig) # hide

fig = Figure() # hide
ax = CairoMakie.Axis(fig[1,1], aspect=1, title = "Initial value of V") # hide
msh = CairoMakie.mesh!(ax, s, color=V, colormap=:jet, colorrange=extrema(V)) # hide
Colorbar(fig[1,2], msh) # hide
save("initial_V.png", fig) # hide

F₁ = map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end

F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  B = 1.0, # Boundary value
  F = t -> t ≥ 1.1 ? F₁ : F₂)
nothing # hide
```

![Initial U Conditions](initial_U.png)
![Initial V Conditions](initial_V.png)

```@example DEC
# Find left and right vertices of mesh
min_x = minimum(x -> x[1], s[:point])
max_x = maximum(x -> x[1], s[:point])
left_wall_idxs = findall(x -> x[1] == min_x, s[:point])
right_wall_idxs = findall(x -> x[1] == max_x, s[:point])
wall_idxs = vcat(left_wall_idxs, right_wall_idxs)

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :rlb => (x,y) -> begin
      x[wall_idxs] .= y
      x
    end
    _ => error("Unmatched operator $my_symbol")
  end
end

fig = Figure() # hide
ax = CairoMakie.Axis(fig[1,1], aspect=1, title = "Highlighted Boundary") # hide
value = zeros(nv(sd)) # hide
value[wall_idxs] .= 1.0 # hide
msh = CairoMakie.mesh!(ax, s, color=value, colormap=:jet, colorrange=(0,2)) # hide
Colorbar(fig[1,2], msh) # hide
save("boundary.png", fig) # hide
nothing # hide
```

![Boundary Visualized](boundary.png)


## Generate the Simulation

We generate our simulation code and store the function in `fₘ` and then run our simulation for `t=11.5` simulated time units.
```@example DEC
sim = evalsim(Brusselator)
fₘ = sim(sd, generate, DiagonalHodge())

u₀ = ComponentArray(U=U, V=V)

tₑ = 11.5

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode
```

## Visualize

We can then use Makie to visualize the evolution of our Brusselator model.
```@example DEC
function save_dynamics(save_file_name)
  time = Observable(0.0)
  u = @lift(soln($time).U)
  v = @lift(soln($time).V)
  f = Figure(; size = (500, 850))
  ax_U = CairoMakie.Axis(f[1,1], title = @lift("Concentration of U at Time $($time)"))
  ax_V = CairoMakie.Axis(f[2,1], title = @lift("Concentration of V at Time $($time)"))

  msh_U = mesh!(ax_U, s, color=u, colormap=:jet, colorrange=extrema(soln(tₑ).U))
  Colorbar(f[1,2], msh_U)

  msh_V = mesh!(ax_V, s, color=v, colormap=:jet, colorrange=extrema(soln(tₑ).V))
  Colorbar(f[2,2], msh_V)

  timestamps = range(0, tₑ, step=1e-1)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("brusselator.gif")
nothing # hide
```

![Brusselator_results_flat](brusselator.gif)

```@example INFO
DocInfo.get_report(info) # hide
```

