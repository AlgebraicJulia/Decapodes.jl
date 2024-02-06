# Introduction to Decapodes

Discrete Exterior Calculus Applied to Partial and Ordinary Differential
Equations (Decapodes) is a diagrammatic language used to express systems of
ordinary and partial differential equations. The Decapode provides a visual
framework for understanding the coupling between variables within a PDE or ODE
system, and a combinatorial data structure for working with them. Below, we
provide a high-level overview of how Decapodes can be generated and interpreted.

## Your First Decapode

We begin with the most basic Decapode, one which only includes a single
variable. In the Decapode graphical paradigm, nodes represent variables and
arrows represent operators which relate variables to each other. Since the
Decapode applies this diagrammatic language specifically to the Discrete
Exterior Calculus (DEC), variables are typed by the dimension and orientation
of the information they contain. So a variable of type `Form0` will be the
0-dimensional data points on some space, or in a discrete context, the
values defined on points of a mesh. Similarly, `Form1` will be values
stored on edges of the mesh, and `Form2` will be values stored on the
surfaces of the mesh. Below, we provide a very simple Decapode with just a
single variable `C`.
and define a convenience function for visualization in later examples.

```@example DEC
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using Catlab, Catlab.Graphics

Variable = @decapode begin
  C::Form0
end;

to_graphviz(Variable)
```

The resulting diagram contains a single node, showing the single variable in
this system. We can then add a second variable:

```@example DEC
TwoVariables = @decapode begin
  C::Form0
  dC::Form1
end;

to_graphviz(TwoVariables)
```

And then can add some relationship between them. In this case, we make an
equation which states that `dC` is the discrete derivative of `C`:

```@example DEC
Equation = @decapode begin
  C::Form0
  dC::Form1

  dC == d₀(C)
end;

to_graphviz(Equation)
```

Here, the two nodes represent the two variables, and the arrow between them
shows how they are related by the discrete derivative.

##  A Little More Complicated

Now that we've seen how to construct a simple equation, it's time to move on to
some actual PDE systems! One classic PDE example is the diffusion equation.
This equation states that the change of concentration at each point is
proportional to the Laplacian of the concentration.

```@example DEC
Diffusion = @decapode begin
  (C, Ċ)::Form0
  ϕ::Form1

  # Fick's first law
  ϕ ==  k(d₀(C))
  # Diffusion equation
  Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
  ∂ₜ(C) == Ċ
end;

to_graphviz(Diffusion)
```

The resulting Decapode shows the relationships between the three variables with
the triangle diagram. Note that these diagrams are automatically layed-out by Graphviz.

## Bring in the Dynamics

Now that we have a reasonably complex PDE, we can demonstrate some of the
developed tooling for actually solving the PDE. Currently, the tooling will
automatically generate an explicit method for solving the system (using
DifferentialEquations.jl to handle time-stepping and instability detection).

We begin this process by importing a mesh. The mesh has been pre-generated
within CombinatorialSpaces, and is generated such that it has periodic boundary
conditions. We will also upload a non-periodic mesh for the sake of
visualization, as well as a mapping between the points on the periodic and
non-periodic meshes.

See CombinatorialSpaces.jl for mesh construction and importing utilities.

```@example DEC
using Catlab.CategoricalAlgebra
using CombinatorialSpaces, CombinatorialSpaces.DiscreteExteriorCalculus
using CairoMakie

plot_mesh = loadmesh(Rectangle_30x10())
periodic_mesh = loadmesh(Torus_30x10())
point_map = loadmesh(Point_Map())

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
wireframe!(ax, plot_mesh)
fig
```

With the mesh uploaded, we also need to convert the Decapode into something
which can be scheduled with explicit time stepping. In order to do this, we take
every variable which is the time derivative of another variable and trace back
the operations needed to compute this. This process essentially generates a computation
graph in the form of a directed wiring diagram.

Since our diagram is already defined, we just need to define a function which implements each of
these symbolic operators and pass them to a scheduler for generating the
function.

Note that we chose to define `k` as a function that multiplies by a value `k`. We could have alternately chosen to represent `k` as a Constant that we multiply by in the Decapode itself.

```@example DEC
using MLStyle

function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :k => x -> 0.05*x
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

Next, we generate the simulation function using `gen_sim` and set up our
initial conditions for this problem.

```@example DEC
sim = eval(gensim(Diffusion))
fₘ = sim(periodic_mesh, generate, DiagonalHodge())

using Distributions
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
mesh!(ax, plot_mesh; color=c[point_map])
fig
```

Finally, we solve this PDE problem using the `Tsit5()` solver and generate an animation of the result!

```@example DEC
using LinearAlgebra
using ComponentArrays
using OrdinaryDiffEq

u₀ = ComponentArray(C=c)

prob = ODEProblem(fₘ, u₀, (0.0, 100.0))
sol = solve(prob, Tsit5());

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t).C[point_map] for t in times]

# Initial frame
fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, plot_mesh; color=colors[1], colorrange = extrema(vcat(colors...)))
Colorbar(fig[1,2], pmsh)
framerate = 30

# Animation
record(fig, "diffusion.gif", range(0.0, 100.0; length=150); framerate = 30) do t
  pmsh.color = sol(t).C[point_map]
end
```

![](diffusion.gif)

## Merging Multiple Physics

Now that we've seen the basic pipeline, it's time for a more complex example
that demonstrates some of the benefits reaped from using Catlab.jl as the
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

  # Fick's first law
  ϕ ==  k(d₀(C))
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
true # hide
```

```@example DEC
to_graphviz(Diffusion)
```

```@example DEC
to_graphviz(Advection)
```

```@example DEC
to_graphviz(Superposition)
```

Next, we define the pattern of composition which we want to compose these
physics under. This pattern of composition is described by an undirected wiring
diagram, which has the individual physics as nodes and the shared variables as
the small junctions.

```@example DEC
compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

to_graphviz(compose_diff_adv, box_labels=:name, junction_labels=:variable, prog="circo")
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

Similar to before, this physics can be compiled and executed. Note that this process
now requires another value to be defined, namely the velocity vector field. We do
this using a custom operator called `flat_op`. This operator is basically the flat
operator from CombinatorialSpaces.jl, but specialized to account for the periodic mesh.

We could instead represent the domain as the surface of an object with equivalent boundaries in 3D.

``` @setup DEC
function closest_point(p1, p2, dims)
    p_res = collect(p2)
    for i in 1:length(dims)
        if dims[i] != Inf
            p = p1[i] - p2[i]
            f, n = modf(p / dims[i])
            p_res[i] += dims[i] * n
            if abs(f) > 0.5
                p_res[i] += sign(f) * dims[i]
            end
        end
    end
    Point3{Float64}(p_res...)
end

function flat_op(s::AbstractDeltaDualComplex2D, X::AbstractVector; dims=[Inf, Inf, Inf])
  tri_map = Dict{Int,Int}(triangle_center(s,t) => t for t in triangles(s))

  map(edges(s)) do e
    p = closest_point(point(s, tgt(s,e)), point(s, src(s,e)), dims)
    e_vec = (point(s, tgt(s,e)) - p) * sign(1,s,e)
    dual_edges = elementary_duals(1,s,e)
    dual_lengths = dual_volume(1, s, dual_edges)
    mapreduce(+, dual_edges, dual_lengths) do dual_e, dual_length
      X_vec = X[tri_map[s[dual_e, :D_∂v0]]]
      dual_length * dot(X_vec, e_vec)
    end / sum(dual_lengths)
  end
end
```

```@example DEC
using LinearAlgebra
using MLStyle

function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :k => x -> 0.05*x
    :∧₀₁ => (x,y) -> begin
      ∧(Tuple{0,1}, sd, x,y)
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

sim = eval(gensim(DiffusionAdvection))
fₘ = sim(periodic_mesh, generate, DiagonalHodge())

velocity(p) = [-0.5, -0.5, 0.0]
v = flat_op(periodic_mesh, DualVectorField(velocity.(periodic_mesh[triangle_center(periodic_mesh),:dual_point])); dims=[30, 10, Inf])

u₀ = ComponentArray(C=c,V=v)

prob = ODEProblem(fₘ, u₀, (0.0, 100.0))
sol = solve(prob, Tsit5());

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t).C[point_map] for t in times]

# Initial frame
fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, plot_mesh; color=colors[1], colorrange = extrema(vcat(colors...)))
Colorbar(fig[1,2], pmsh)
framerate = 30

# Animation
record(fig, "diff_adv.gif", range(0.0, 100.0; length=150); framerate = 30) do t
  pmsh.color = sol(t).C[point_map]
end
```

![](diff_adv.gif)
```
