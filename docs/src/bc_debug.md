# Simulation Setup

This tutorial showcases some of the other features included in the Decapodes.jl
package. Currently, these features are the treatment of boundary conditions and
the simulation debugger interface. To begin, we set up the same
advection-diffusion problem presented in the Overview section.
As before, we define the Diffusion, Advection, and Superposition components,
and now include a BC (Bounday Condition) component. Decapodes.jl interprets any
`Hom` which begins with a `∂` as a boundary condition. These boundary
conditions recieve special treatment at the scheduling step. Below we show the
graphical rendering of this boundary condition diagram, which we will use to
impose a Dirichlet condition on the time derivative of concentration at the
mesh boundary.

```@example Debug
using Catlab
using Catlab.Graphics
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes

Diffusion = @decapode begin
  C::Form0
  ϕ::Form1

  # Fick's first law
  ϕ ==  k(d₀(C))
end

Advection = @decapode begin
  C::Form0
  ϕ::Form1
  V::Constant

  ϕ == ∧₀₁(C,V)
end

Superposition = @decapode begin
  (C, C_up)::Form0
  (ϕ, ϕ₁, ϕ₂)::Form1

  ϕ == ϕ₁ + ϕ₂
  C_up == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
end

BoundaryConditions = @decapode begin
  (C, C_up)::Form0

  # Temporal boundary
  ∂ₜ(C) == Ċ

  # Spatial boundary
  Ċ == ∂C(C_up)
end

to_graphviz(BoundaryConditions)
```

As before, we compose these physics components over our wiring diagram.


```@example Debug
compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  bc(C, C_up)
  superposition(ϕ₁, ϕ₂, ϕ, C_up, C)
end

to_graphviz(compose_diff_adv, box_labels=:name, junction_labels=:variable, prog="circo")
```

```@example Debug
DiffusionAdvection_cospan = oapply(compose_diff_adv,
                  [Open(Diffusion, [:C, :ϕ]),
                   Open(Advection, [:C, :ϕ, :V]),
                   Open(BoundaryConditions, [:C, :C_up]),
                   Open(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C_up, :C])])
DiffusionAdvection = apex(DiffusionAdvection_cospan)

to_graphviz(DiffusionAdvection)
```

When this is scheduled, Decapodes will apply any boundary conditions
immediately after the impacted value is computed. This implementation choice
ensures that this boundary condition holds true for any variables dependent on
this variable, though also means that the boundary conditions on a variable
have no immediate impact on the variables this variable is dependent on.

In the visualization below, wee see that the final operation
executed on the data is the boundary condition we are enforcing on the change
in concentration.


```@example Debug
to_graphviz(DiffusionAdvection)
```

Next we import the mesh we will use. In this case, we are wanting to impose
boundary conditions and so we will use the `plot_mesh` from the previous
example instead of the mesh with periodic boundary conditions. Because the mesh
is only a primal mesh, we also generate and subdivide the dual mesh.

`Rectangle_30x10` is a default mesh that is downloaded via `Artifacts.jl` when a user installs Decapodes. Via CombinatorialSpaces.jl, we can instantiate any `.obj` file of triangulated faces as a simplicial set.

```@example Debug
using CombinatorialSpaces, CombinatorialSpaces.DiscreteExteriorCalculus
using CairoMakie

plot_mesh = loadmesh(Rectangle_30x10())

# Generate the dual mesh
plot_mesh_dual = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}}(plot_mesh)
# Calculate distances and subdivisions for the dual mesh
subdivide_duals!(plot_mesh_dual, Circumcenter())

fig = Figure()
ax = Axis(fig[1,1], aspect = AxisAspect(3.0))
wireframe!(ax, plot_mesh)
fig
```

Finally, we define our operators, generate the simulation function, and compute
the simulation. Note that when we define the boundary condition operator, we
hardcode the boundary indices and values into the operator itself. We also move
the initial concentration to the left, so that we are able to see a constant
concentration on the left boundary which will act as a source in the flow. The
modified initial condition is shown below:

```@example Debug
using LinearAlgebra
using ComponentArrays
using MLStyle
using CombinatorialSpaces.DiscreteExteriorCalculus: ∧
include("../../examples/boundary_helpers.jl")

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :k => x -> 0.05*x
    :∂C => x -> begin
      boundary = boundary_inds(Val{0}, sd)
      x[boundary] .= 0
      x
    end
    :∧₀₁ => (x,y) -> begin
      ∧(Tuple{0,1}, sd, x,y)
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

using Distributions
c_dist = MvNormal([1, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in plot_mesh_dual[:point]]

fig = Figure()
ax = Axis(fig[1,1], aspect = AxisAspect(3.0))
mesh!(ax, plot_mesh; color=c[1:nv(plot_mesh)])
fig
```

And the simulation result is then computed and visualized below.

```@example Debug
using OrdinaryDiffEq

sim = eval(gensim(DiffusionAdvection))
fₘ = sim(plot_mesh_dual, generate)

velocity(p) = [-0.5, 0.0, 0.0]
v = ♭(plot_mesh_dual, DualVectorField(velocity.(plot_mesh_dual[triangle_center(plot_mesh_dual),:dual_point]))).data

u₀ = ComponentArray(C=c)
params = (V = v,)

prob = ODEProblem(fₘ, u₀, (0.0, 100.0), params)
sol = solve(prob, Tsit5());

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t).C for t in times]
extrema
# Initial frame
fig = Figure()
ax = Axis(fig[1,1], aspect = AxisAspect(3.0))
pmsh = mesh!(ax, plot_mesh; color=colors[1], colorrange = extrema(vcat(colors...)))
Colorbar(fig[1,2], pmsh)
framerate = 30

# Animation
record(fig, "diff_adv_right.gif", range(0.0, 100.0; length=150); framerate = 30) do t
  pmsh.color = sol(t).C
end
```

![](diff_adv_right.gif)
