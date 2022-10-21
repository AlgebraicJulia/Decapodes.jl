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
using Decapodes, Decapodes.Diagrams, Decapodes.OpenDiagrams
using Catlab.Present, Catlab.Graphics, Catlab.Programs

@present DiffusionQuantities <: Decapodes2D begin
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  ∂C::Hom(Form0(X), Form0(X))   # concentration boundary condition
end;

Diffusion = @decapode DiffusionQuantities begin
  C::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
end

Advection = @decapode DiffusionQuantities begin
  C::Form0{X}
  (V, ϕ)::Form1{X}
  ϕ == ∧₀₁{X}(C,V)
end

Superposition = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  (ϕ, ϕ₁, ϕ₂)::Form1{X}

  ϕ == ϕ₁ + ϕ₂
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end

BoundaryConditions = @decapode DiffusionQuantities begin
  (Ċ, Ċ_bound)::Form0{X}
  ∂C(Ċ) == Ċ_bound
end

draw_diagram(BoundaryConditions)
```

As before, we compose these physics components over our wiring diagram.


```@example Debug
compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  bc(Ċ)
  superposition(ϕ₁, ϕ₂, ϕ, Ċ, C)
end

to_graphviz(compose_diff_adv, box_labels=:name, junction_labels=:variable,
            graph_attrs=Dict(:start => "2"))
```



```@example Debug
DiffusionAdvection = oapply(compose_diff_adv,
                  [OpenDiagram(Diffusion, [:C, :ϕ]),
                   OpenDiagram(Advection, [:C, :ϕ, :V]),
                   OpenDiagram(BoundaryConditions, [:Ċ]),
                   OpenDiagram(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :Ċ, :C])])

draw_diagram(DiffusionAdvection.functor)
```

When this is scheduled, Decapodes will apply any boundary conditions
immediately after the impacted value is computed. This implementation choice
ensures that this boundary condition holds true for any variables dependent on
this variable, though also means that the boundary conditions on a variable
have no immediate impact on the variables this variable is dependent on.

Below, we see the generated schedule, which shows that the final operation
executed on the data is the boundary condition we are enforcing on the change
in concentration.


```@example Debug
using Decapodes.Schedules

explicit_ts = diag2dwd(DiffusionAdvection.functor)
to_graphviz(explicit_ts, orientation=LeftToRight)
```

Next we import the mesh we will use. In this case, we are wanting to impose
boundary conditions and so we will use the `plot_mesh` from the previous
example instead of the mesh with periodic boundary conditions. Because the mesh
is only a primal mesh, we also generate and subdivide the dual mesh.


```@example Debug
using Catlab.CategoricalAlgebra
using CombinatorialSpaces, CombinatorialSpaces.DiscreteExteriorCalculus
using CairoMakie

plot_mesh = loadmesh(Rectangle_30x10())

# Generate the dual mesh
plot_mesh_dual = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}}(plot_mesh)
# Calculate distances and subdivisions for the dual mesh
subdivide_duals!(plot_mesh_dual, Circumcenter())


fig, ax, ob = wireframe(plot_mesh)
ax.aspect = AxisAspect(3.0)
fig
```

Finally, we define our operators, generate the simulation function, and compute
the simulation. Note that when we define the boudary condition operator, we
hardcode the boundary indices and values into the operator itself. We also move
the initial concentration to the left, so that we are able to see a constant
concentration on the left boundary which will act as a source in the flow. The
modified initial condition is shown below:


```@example Debug
using LinearAlgebra
using Decapodes.Examples, Decapodes.Simulations
using CombinatorialSpaces.DiscreteExteriorCalculus: ∧

funcs = sym2func(plot_mesh_dual)

funcs[:k] = Dict(:operator => 0.05 * I(ne(plot_mesh_dual)), :type => MatrixFunc())
funcs[:⋆₁] = Dict(:operator => ⋆(Val{1}, plot_mesh_dual, hodge=DiagonalHodge()),
                  :type => MatrixFunc());
funcs[:∧₀₁] = Dict(:operator => (r, c,v)->r .= ∧(Tuple{0,1}, plot_mesh_dual, c, v), :type => InPlaceFunc())

boundary = Examples.boundary_inds(Val{0}, plot_mesh)
funcs[:∂C] = Dict(:operator => (∂ċ, ċ)->(∂ċ .= ċ; ∂ċ[boundary] .= 0), :type => InPlaceFunc())


using Distributions
c_dist = MvNormal([1, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in plot_mesh_dual[:point]]

fig, ax, ob = mesh(plot_mesh; color=c[1:nv(plot_mesh)])
ax.aspect = AxisAspect(3.0)
fig
```

And the simulation result is then computed and visualized below.


```@example Debug
using OrdinaryDiffEq
func, code = gen_sim(explicit_ts, funcs, plot_mesh_dual; autodiff=false, params = [:V]);

velocity(p) = [-0.5, 0.0, 0.0]
v = ♭(plot_mesh_dual, DualVectorField(velocity.(plot_mesh_dual[triangle_center(plot_mesh_dual),:dual_point]))).data

prob = ODEProblem(func, c, (0.0, 100.0))
sol = solve(prob, Tsit5(), p=v);

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t)[1:nv(plot_mesh)] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "diff_adv_right.gif", range(0.0, 100.0; length=150); framerate = 30) do t
ob.color = sol(t)[1:nv(plot_mesh)]
end
```

![](diff_adv_right.gif)

# Debug the Simulation

The benefit of having the computation graph which we generated above means that
we have the opportunity to inspect the simulation at different stages of
computation. Since each wire in the computation diagram has the data of a
particular form on the mesh, this data can be visualized. The first step of
getting at this data, though, is understanding the index associated with each
wire. This key between indices and wires can be generated with the `sim_key`
function.


```@example Debug
using Decapodes.Debug

sim_key(explicit_ts)
```

Another key piece of data is what the initial state is at the time-step you are
debugging. We choose the time step `t = 10` for this example, shown below.


```@example Debug
t = 10.0
fig, ax, ob = mesh(plot_mesh; color=sol(t)[1:nv(plot_mesh_dual)])
xlims!(ax, (0,10.0))
fig
```

Now we can see what the value is of the result of the product between velocity
and the concentration.


```@example Debug
fig, ax, ob = draw_wire(plot_mesh, plot_mesh_dual, explicit_ts, func, sol(10.0), 2;p = v, xlim=(0, 10.0), ylim=(0.0,10.0))
fig
```

Now, this information doesn't seem all that useful at first glance. The result
is basically a heatmap of the magnitude of the vector field across each edge of
the mesh. Here, we see that this product has a greater magnitude on edges that
are both aligned with the flow and have a high concentraiton at that point.

Something which is a little more useful is to includ arrows! By using the
`use_arrows` keyword argument, we are able to get arrows, directed along the
mesh edge elements, which show the direction of the vectorfield at that point.
The color of these arrows shows the magnitude of the vectorfield in that
direction. Note that this method of plotting arrows with CairoMakie is fairly
computationally expensive, so as the `n_arrows` argument increases, so will the
time it takes to render.


```@example Debug
fig, ax, ob = draw_wire(plot_mesh, plot_mesh_dual, explicit_ts, func, sol(10.0), 2;p = v, xlim=(0, 10.0), ylim=(0.0,10.0), use_arrows=true, n_arrows=1200)
fig
```

Note that all of the arrows from the result of the product are pointing in the
direction of flow. This intuitively makes sense for how the concentration flow
vectorfield should look.

Below, we show the concentration flow resulting just from the gradient.


```@example Debug
fig, ax, ob = draw_wire(plot_mesh, plot_mesh_dual, explicit_ts, func, sol(10.0), 7;p = v, xlim=(0, 10.0), ylim=(0.0,10.0), use_arrows=true, n_arrows=1200)
fig
```

Here, we see that the arrows are pointing away from areas of high concentration
to areas of low concentration, and that the magnitude of the arrows grows
greater as the slope of concentration is greater.

Next, we can show what happens after the effects of diffusion and advection are
added together.


```@example Debug
fig, ax, ob = draw_wire(plot_mesh, plot_mesh_dual, explicit_ts, func, sol(10.0), 4;p = v, xlim=(0, 10.0), ylim=(0.0,10.0), use_arrows=true, n_arrows=1200)
fig
```

Finally, in order to compute the change of concentration, we visualize what the
resulting change in concentration is.


```@example Debug

fig, ax, ob = draw_wire(plot_mesh, plot_mesh_dual, explicit_ts, func, sol(10.0), 3;p = v, xlim=(0, 10.0), ylim=(0.0,10.0), use_arrows=true, n_arrows=1200)
fig
```

Here, we can see the positive change which is weighted heavily to the right of
the existing distribution, and the slight negative change which follows from
where the distribution is flowing.

