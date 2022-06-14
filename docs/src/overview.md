# Introduction to DECAPODEs

Discrete Exterior Calculus Applied to Partial and Ordinary Differential
Equations (DECAPODE) is a diagrammatic language used to express systems of
ordinary and partial differential equations. The DECAPODE provides a visual
framework for understanding the coupling between variables within a PDE or ODE
system, and a combinatorial data structure for working with them. Below, we
provide a high-level overview of how DECAPODEs can be generated and interpreted.

## Your First DECAPODE

We begin with the most basic DECAPODE, one which only includes a single
variable. In the DECAPODE graphical paradigm, nodes represent variables and
arrows represent operators which relate variables to each other. Since the
DECAPODE applies this diagrammatic language specifically to the Discrete
Exterior Calculus (DEC), variables are typed by the dimension and orientation
of the information they contain. So a variable of type `Form0{X}` will be the
0-dimensional data points on the space `X`, or in a discrete context, the
values defined on points of a mesh `X`. Similarly, `Form1{X}` will be values
stored on edges of the mesh, and `Form2{X}` will be values stored on the
surfaces of the mesh. Below, we provide a very simple DECAPODE with just a
single variable `C`. In this example, we also provide the necessary imports,
and define a convenience function for visualization in later examples.
```@example DEC
using Decapodes, Decapodes.Diagrams
using Catlab.Present, Catlab.Graphics

Variable = @decapode Decapodes2D begin
  C::Form0{X}
end;

draw_equation(decapode) = to_graphviz(decapode, node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"4.0"))

draw_equation(Variable)
```

The resulting diagram contains a single node, showing the single variable in
this system. We can then add a second variable:

```@example DEC
TwoVariables = @decapode Decapodes2D begin
  C::Form0{X}
  dC::Form1{X}
end;

draw_equation(TwoVariables)
```

And then can add some relationship between them. In this case, we make an
equation which states that `dC` is the discrete derivative of `C`:

```@example DEC
Equation = @decapode Decapodes2D begin
  C::Form0{X}
  dC::Form1{X}

  dC == d₀{X}(C)
end;

draw_equation(Equation)
```

Here, the two nodes represent the two variables, and the arrow between them
shows how they are related by the discrete derivative.

##  A Little More Complicated

Now that we've seen how to construct a simple equation, it's time to move on to
some actual PDE systems! One classic PDE example is the diffusion equation.
This equation states that the change of concentration at each point is
proportional to the laplacian of the concentration. One issue that we run into
here, though, is that there isn't a "proportionality" operator in the default
DECAPODEs syntax `Decapodes2D`. Thus, in this next example, we will first
extend the `Decapodes2D` syntax and then define the DECAPODE for diffusion.

```@example DEC
@present DiffusionQuantities <: Decapodes2D begin
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
end;

Diffusion = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
  # Diffusion equation
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end;

draw_equation(Diffusion)
```

The resulting DECAPODE shows the relationships between the three variables with
the triangle diagram. Note that automatic layout of these labels can result in
confusion as to which edge each label corresponds, but any confusion can be
resolved by referring back to the original `@decapode` definition.

## Bring in the Dynamics

Now that we have a reasonably complex PDE, we can demonstrate some of the
developed tooling for actually solving the PDE. Currently, the tooling will
automatically generate an explicit method for solving the system (using
DifferentialEquatons.jl for the actual computation).

We begin this process by importing a mesh. The mesh has been pre-generated
within CombinatorialSpaces, and is generated such that it has periodic boundary
conditions. We will also upload a non-periodic mesh for the sake of
visualization, as well as a mapping between the points on the periodic and
non-periodic meshes.
```@example DEC
using Catlab.CategoricalAlgebra
using CombinatorialSpaces, CombinatorialSpaces.DiscreteExteriorCalculus
using CairoMakie
using JSON
using HTTP: get

@info pwd() #hide
periodic_mesh = parse_json_acset(EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}},
                                 String(get("https://raw.githubusercontent.com/AlgebraicJulia/Decapodes.jl/main/docs/assets/meshes/periodic_mesh.json").body))
plot_mesh = parse_json_acset(EmbeddedDeltaSet2D{Bool, Point3{Float64}},
                                 String(get("https://raw.githubusercontent.com/AlgebraicJulia/Decapodes.jl/main/docs/assets/meshes/plot_mesh.json").body))
point_map = JSON.parse(String(get("https://raw.githubusercontent.com/AlgebraicJulia/Decapodes.jl/main/docs/assets/meshes/point_map.json").body))

fig, ax, ob = wireframe(plot_mesh)
ax.aspect = AxisAspect(3.0)
fig
```

With the mesh uploaded, we also need to convert the DECAPODE into something
which can be scheduled as an explicit time step. In order to do this, we take
every variable which is the time derivative of another variable and trace back
the operations needed to compute this. This process generates a computation
graph in the form of a directed wiring diagram, as shown below.

```@example DEC
using Decapodes.Schedules

explicit_ts = diag2dwd(Diffusion)
to_graphviz(explicit_ts, orientation=LeftToRight)
```

As can be seen, diffusion has a very simple explicit time step. With this
diagram defined, we just need to define a function which implements each of
these symbolic operators and pass them to a scheduler for generating the
function. The basic DEC operators can be computed with the `sym2func` operator,
and more complex or custom functions (like the proportionality constant `k`)
can be defined manually. We choose to also define the hodge star (⋆₁) on primal
1-forms manually, as this mesh allows for a diagonal hodge star which is
significantly more efficient for computation.

```@example DEC
using LinearAlgebra
using Decapodes.Examples, Decapodes.Simulations

funcs = sym2func(periodic_mesh)

funcs[:k] = Dict(:operator => 0.05 * I(ne(periodic_mesh)), :type => MatrixFunc())
funcs[:⋆₁] = Dict(:operator => ⋆(Val{1}, periodic_mesh, hodge=DiagonalHodge()),
                  :type => MatrixFunc());
```

Next, we generate the simulation function using `gen_sim` and set up our
initial conditions for this problem.

```@example DEC
func, code = gen_sim(explicit_ts, funcs, periodic_mesh; autodiff=false);

using Distributions
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

fig, ax, ob = mesh(plot_mesh; color=c[point_map])
ax.aspect = AxisAspect(3.0)
fig
```

Finally, we solve this PDE problem using the `Tsit5()` solver and generate an animation of the result!

```@example DEC
using DifferentialEquations

prob = ODEProblem(func, c, (0.0, 100.0))
sol = solve(prob, Tsit5());

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t)[point_map] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "diffusion.gif", range(0.0, 100.0; length=150); framerate = 30) do t
ob.color = sol(t)[point_map]
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
```
```@example DEC
draw_equation(Diffusion)
```
```@example DEC
draw_equation(Advection)
```
```@example DEC
draw_equation(Superposition)
```

Next, we define the pattern of composition which we want to compose these
physics under. This pattern of composition is described by an undirected wiring
diagram, which has the individual physics as nodes and the shared variables as
the small junctions.

```@example DEC
using Catlab.Programs

compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

to_graphviz(compose_diff_adv, box_labels=:name, junction_labels=:variable,
            graph_attrs=Dict(:start => "2"))

```

After this, the physics can be composed as follows:

```@example DEC
using Decapodes.OpenDiagrams
DiffusionAdvection = oapply(compose_diff_adv,
                  [OpenDiagram(Diffusion, [:C, :ϕ]),
                   OpenDiagram(Advection, [:C, :ϕ, :V]),
                   OpenDiagram(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])])

draw_equation(DiffusionAdvection.functor)
```

Similar to before, this physics can be converted to a directed wiring diagram,
compiled, and executed. Note that this process now requires another value to be
defined, namely the velocity vector field. We do this using a custom operator
called `flat_op`. This operator is basically the flat operator from
CombinatorialSpaces.jl, but specialized to account for the periodic mesh.

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
  # XXX: Creating this lookup table shouldn't be necessary. Of course, we could
  # index `tri_center` but that shouldn't be necessary either. Rather, we should
  # loop over incident triangles instead of the elementary duals, which just
  # happens to be inconvenient.
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
explicit_ts = diag2dwd(DiffusionAdvection.functor)
to_graphviz(explicit_ts, orientation=LeftToRight)
```

```@example DEC
using CombinatorialSpaces.DiscreteExteriorCalculus: ∧
funcs[:∧₀₁] = Dict(:operator => (r, c,v)->r .= ∧(Tuple{0,1}, periodic_mesh, c, v), :type => InPlaceFunc())

func, code = gen_sim(explicit_ts, funcs, periodic_mesh; autodiff=false, params = [:V]);

velocity(p) = [-0.5, -0.5, 0.0]
v = flat_op(periodic_mesh, DualVectorField(velocity.(periodic_mesh[triangle_center(periodic_mesh),:dual_point])); dims=[30, 10, Inf])

prob = ODEProblem(func, c, (0.0, 100.0))
sol = solve(prob, Tsit5(), p=v);

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t)[point_map] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "diff_adv.gif", range(0.0, 100.0; length=150); framerate = 30) do t
ob.color = sol(t)[point_map]
end
```

![](diff_adv.gif)
```
