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
end

to_graphviz(Variable)
```

The resulting diagram contains a single node, showing the single variable in
this system. The arrow into this node represents that `C` is a source term.

Note that to display these diagrams you will require [Graphviz](https://graphviz.org/).

We can add a second variable.

```@example DEC
TwoVariables = @decapode begin
  C::Form0
  dC::Form1
end

to_graphviz(TwoVariables)
```

We can also add a relationship between them. In this case, we make an
equation which states that `dC` is the derivative of `C`.

```@example DEC
Equation = @decapode begin
  C::Form0
  dC::Form1

  dC == d(C)
end

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

end

to_graphviz(Diffusion)
```

The resulting Decapode shows the relationships between the three variables with
the triangle diagram.

## Bring in the Dynamics

Now that we have a reasonably complex PDE, we can demonstrate some of the
developed tooling for actually solving the PDE. Currently, the tooling will
automatically generate an explicit method for solving the system (using
[OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) to handle time-stepping and instability detection).

To start with, we solve our PDE on a flat `triangulated_grid` as provided by [CombinatorialSpaces.jl](https://algebraicjulia.github.io/CombinatorialSpaces.jl/stable/meshes/#CombinatorialSpaces.Meshes.triangulated_grid).

```@example DEC
using CairoMakie
using CombinatorialSpaces
using GeometryBasics: Point3

mesh = triangulated_grid(30,10,1,1,Point3{Float64})

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
wireframe!(ax, mesh)
fig
```

We first create a mesh that is "dual" to the original mesh by creating an `EmbeddedDeltaDualComplex2D` and calling `subdivide_duals!`. This is done to allow important DEC operators to be created and used in the simulation.

We then compile the simulation by using `gensim` and create functional simulation by calling the evaluated `sim`  with the mesh. Decapodes supports providing custom functions to the simulator, but for the time being we will keep things simple by passing in `nothing` to this value.

```@example DEC
dualmesh = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}}(mesh)
subdivide_duals!(dualmesh, Circumcenter())
sim = eval(gensim(Diffusion))
fₘ = sim(dualmesh, nothing, DiagonalHodge())
```

We go ahead and set up our initial conditions for this problem. In this case we generate a Gaussian over our mesh.

```@example DEC
using Distributions
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in mesh[:point]]

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
mesh!(ax, mesh; color=c)
fig
```

Finally, we solve this PDE problem using the `Tsit5()` solver provided by [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/explicit/Tsit5/). We will need to make two ComponentArrays to store the state variable data and constants, respectively.

```@example DEC
using LinearAlgebra
using ComponentArrays
using OrdinaryDiffEq

u₀ = ComponentArray(C=c)
constants = ComponentArray(k=0.05)

tₑ = 100
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), constants)
soln = solve(prob, Tsit5())
soln.retcode
```

Now that the simulation has succeeded we can plot out our results with [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl).

```@example DEC
fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, mesh; color=c, colorrange = extrema(c))
Colorbar(fig[1,2], pmsh)


times = range(0.0, tₑ, length=150)
record(fig, "diffusion.gif", times; framerate = 30) do t
  pmsh.color = soln(t).C
end
nothing # hide
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

For this example, we will invoke `Catlab`'s exported [`@relation`](https://algebraicjulia.github.io/Catlab.jl/v0.9/apis/programs/#Catlab.Programs.RelationalPrograms.@relation-Tuple) to specify a composition pattern. We are essentially saying here that the `C` in each Decapode should be fused into the same variable. Likewise for `ϕ₁` between `diffusion` and `superposition` as well as `ϕ₂` between `advection` and `superposition`.

```@example DEC
using Catlab
compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

draw_composition(compose_diff_adv)
```

After this, the physics can be composed as follows. Here we specify, for example, that the `ϕ₁` that we wish to have in our final Decapode (as declared above) is called `ϕ` in `diffusion` and `ϕ₁` in `superposition`. Notice the fused result in `DiffusionAdvection`.

```@example DEC
DiffusionAdvection_cospan = oapply(compose_diff_adv,
                  [Open(Diffusion, [:C, :ϕ]),
                   Open(Advection, [:C, :ϕ, :V]),
                   Open(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])])
DiffusionAdvection = apex(DiffusionAdvection_cospan)

to_graphviz(DiffusionAdvection)
```

Similar to before, this physics can be compiled, executed, and plotted.

Note that for our purposes we require a velocity vector field for the advection (remember we called this term `V` in our Decapode). We define the velocity field so as to move our fluid to the right (the negative is to account for the orientation of the edges of the mesh) and call `eval_constant_primal_form` to convert this field into a `Form1` usable by the Decapode.

We then run the simulation and create video out of the solution.

```@example DEC
using LinearAlgebra
using MLStyle
using StaticArrays
import CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form

sim = eval(gensim(DiffusionAdvection))
fₘ = sim(dualmesh, nothing, DiagonalHodge())

v = eval_constant_primal_form(dualmesh, SVector{3}(-0.25, 0.0, 0.0))

u₀ = ComponentArray(C=c,V=v)
constants = ComponentArray(diffusion_k=0.05)

tₑ = 50
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), constants)
sol = solve(prob, Tsit5())

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, mesh; color=c, colorrange = extrema(c))
Colorbar(fig[1,2], pmsh)

times = range(0.0, tₑ, length=150)
record(fig, "diff_adv.gif", times; framerate = 30) do t
  pmsh.color = soln(t).C
end
nothing # hide
```

![Your first composed Decapode!](diff_adv.gif)

```@example INFO
DocInfo.get_report(info) # hide
```
