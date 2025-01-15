# Introduction to Decapodes

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..", "docinfo.jl"))
info = DocInfo.Info()
```

Discrete Exterior Calculus Applied to Partial and Ordinary Differential
Equations (Decapodes) is a diagrammatic language used to express systems of
ordinary and partial differential equations. Decapodes provides a visual
framework for understanding the coupling between variables within a PDE or ODE
system, and a combinatorial data structure for working with them. Below, we
provide a high-level overview of how Decapodes can be generated and interpreted.

## Your First Decapode

In the Decapodes graphical paradigm, nodes represent variables and
arrows represent operators which relate variables to each other. Since
Decapodes applies this diagrammatic language specifically to the Discrete
Exterior Calculus (DEC), variables are typed by the dimension and orientation
of the information they contain. So a variable of type `Form0` will be the
0-dimensional data points defined the vertices of a mesh. Similarly, `Form1` will be values
stored on edges of the mesh and `Form2` will be values stored on the
surfaces of the mesh.

Below, we provide a Decapode with just a single variable `C` and display it.

```@example DEC
using Decapodes
using DiagrammaticEquations
using Catlab

Variable = @decapode begin
  C::Form0
end;

to_graphviz(Variable)
```

The resulting diagram contains a single node, showing the single variable in
this system. We can add a second variable:

```@example DEC
TwoVariables = @decapode begin
  C::Form0
  dC::Form1
end;

to_graphviz(TwoVariables)
```

We can also add a relationship between them. In this case, we make an
equation which states that `dC` is the derivative of `C`:

```@example DEC
Equation = @decapode begin
  C::Form0
  dC::Form1

  dC == d(C)
end;

to_graphviz(Equation)
```

Here, the two nodes represent the two variables, and the arrow between them
shows how they are related by the derivative.

## A Little More Complicated

Now that we've seen how to construct a simple equation, it's time to move on to
some actual PDE systems! One classic PDE example is the diffusion equation.
This equation states that the change of concentration at each point is
proportional to the Laplacian of the concentration.

```@example DEC
Diffusion = @decapode begin
  (C, Ċ)::Form0
  ϕ::Form1
  k::Constant

  # Fick's first law
  ϕ ==  k*(d₀(C))

  # Diffusion equation
  Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
  ∂ₜ(C) == Ċ

end;

to_graphviz(Diffusion)
```

The resulting Decapode shows the relationships between the three variables with
the triangle diagram. Note that these diagrams are automatically layed-out by [Graphviz](https://graphviz.org/).

## Bring in the Dynamics

Now that we have a reasonably complex PDE, we can demonstrate some of the
developed tooling for actually solving the PDE. Currently, the tooling will
automatically generate an explicit method for solving the system (using
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl?tab=readme-ov-file) to handle time-stepping and instability detection).

`Torus_30x10` is a default mesh that is downloaded via `Artifacts.jl` when a user installs [CombinatorialSpaces.jl](https://github.com/AlgebraicJulia/CombinatorialSpaces.jl). If we wanted, we could also instantiate any `.obj` file of triangulated faces as a simplicial set although we do not here.

We will also upload a non-periodic mesh for the sake of
visualization, as well as a mapping between the points on the periodic and
non-periodic meshes.

```@example DEC
using CairoMakie
using CombinatorialSpaces
using GeometryBasics

mesh = triangulated_grid(20,10,0.1,0.1,Point3{Float64})

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
wireframe!(ax, mesh)
fig
```
We then compile the simulation by using `gensim` and create functional simulation by calling the evaluated `sim`  with the mesh. Decapodes supports providing custom functions to the simulator, but for the time being we will keep things simple by passing in `nothing` to this value. Necessarily, we will create a "dual mesh," an explanation can be found [HERE].

[TODO: subdivide_duals! warrants explanation]

```@example DEC
dualmesh = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}}(mesh)
subdivide_duals!(dualmesh, Barycenter())
sim = eval(gensim(Diffusion))
fₘ = sim(dualmesh, nothing, DiagonalHodge())
```

We go ahead and set up our initial conditions for this problem. In this case we generate a Gaussian and apply it to our mesh.

```@example DEC
using Distributions
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in mesh[:point]]

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
mesh!(ax, mesh; color=c)
fig
```

Finally, we solve this PDE problem using the `Tsit5()` solver provided by DifferentialEquations.jl. We will need to make two ComponentArrays to store the state variable data and constants, respectively.

```@example DEC
using LinearAlgebra
using ComponentArrays
using OrdinaryDiffEq

u₀ = ComponentArray(C=c)
constants = ComponentArray(k=0.5)

prob = ODEProblem(fₘ, u₀, (0.0, 100.0), constants)
soln = solve(prob, Tsit5());
soln.retcode
```

Now that the simulation has succeeded we can plot out our results with [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl).

```@example DEC
# Plot the result
times = range(0.0, 100.0, length=150)

# Initial frame
fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, mesh; color=c, colorrange = extrema(c))
Colorbar(fig[1,2], pmsh)

# Animation
record(fig, "diffusion.gif", times; framerate = 30) do t
  pmsh.color = soln(t).C
end
```

![Your first Decapode!](diffusion.gif)

## Merging Multiple Physics

Now that we've seen the basic pipeline, it's time for a more complex example
that demonstrates some of the benefits reaped from using [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl) as the
backend to our data structures. In this example, we will take two separate
physics (diffusion and advection), and combine them together using a
higher-level composition pattern.

We begin by defining the three systems we need. The first two systems
are the relationships between concentration and flux under diffusion and
advection respectively. The third is the relationship between the two fluxes
and the change of concentration under superposition of fluxes.

```@example DEC
Diffusion = @decapode begin
  C::Form0
  ϕ::Form1
  k::Constant

  # Fick's first law
  ϕ == k*(d₀(C))
end

Advection = @decapode begin
  C::Form0
  ϕ::Form1
  V::Form1

  ϕ == ∧₀₁(C,V)
end

Superposition = @decapode begin
  (C, Ċ)::Form0
  (ϕ, ϕ₁, ϕ₂)::Form1

  ϕ == ϕ₁ + ϕ₂
  Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
  ∂ₜ(C) == Ċ
end
nothing # hide
```

The diffusion Decapode.

```@example DEC
to_graphviz(Diffusion) # hide
```

The advection Decapode.

```@example DEC
to_graphviz(Advection) # hide
```

And the superposition Decapode.

```@example DEC
to_graphviz(Superposition) # hide
```

Next, we define the pattern of composition which we want to compose these
physics under. This pattern of composition is described by an undirected wiring
diagram, which has the individual physics as nodes and the shared variables as
the small junctions.

For this example, we will now invoke `Catlab`, which exports
[`@relation`](https://algebraicjulia.github.io/Catlab.jl/v0.9/apis/programs/#Catlab.Programs.RelationalPrograms.@relation-Tuple)
necessary for specifying a composition pattern.

```@example DEC
using Catlab
compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

draw_composition(compose_diff_adv)
```

After this, the physics can be composed as follows:

```@example DEC
DiffusionAdvection_cospan = oapply(compose_diff_adv,
                  [Open(Diffusion, [:C, :ϕ]),
                   Open(Advection, [:C, :ϕ, :V]),
                   Open(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])])
DiffusionAdvection = apex(DiffusionAdvection_cospan)

to_graphviz(DiffusionAdvection)
```

We used `compose_diff_adv` to specify how physics should be
composed, and then provided a list of physical models (Decapodes) along with
their variables which will be associated to those in the pattern.

Similar to before, this physics can be compiled, executed, and plotted. Note that this process
now requires another value to be defined, namely the velocity vector field. 

TODO: this now slams into the wall. needs fixing.

```@example DEC
using LinearAlgebra
using MLStyle

sim = eval(gensim(DiffusionAdvection))
fₘ = sim(dualmesh, nothing, DiagonalHodge())

velocity(p) = [-0.5, -0.5, 0.0]
v = flat(dualmesh, velocity.(dualmesh[triangle_center(dualmesh),:dual_point]),
DPPFlat())

u₀ = ComponentArray(C=c,V=v)
constants = ComponentArray(diffusion_k=0.5)

prob = ODEProblem(fₘ, u₀, (0.0, 100.0), constants)
sol = solve(prob, Tsit5());

# Plot the result

# Initial frame
fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, mesh; color=c, colorrange = extrema(c))
Colorbar(fig[1,2], pmsh)

# Animation
times = range(0.0, 100.0, length=150)
record(fig, "diff_adv.gif", times; framerate = 30) do t
  pmsh.color = sol(t).C
end
```

![Your first composed Decapode!](diff_adv.gif)

```@example INFO
DocInfo.get_report(info) # hide
```
