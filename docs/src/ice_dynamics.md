# Halfar's model of glacial flow

Let's model glacial flow using a model of how ice height of a glacial sheet changes over time, from P. Halfar's 1981 paper: "On the dynamics of the ice sheets".

```
# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes

# External Dependencies
using MLStyle
using MultiScaleArrays
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using SparseArrays
using Statistics
using GLMakie
using GeometryBasics: Point2
Point2D = Point2{Float64};
```

# Defining the models

The first step is to find a suitable equation for our model, and translate it into the Discrete Exterior Calculus. The Exterior Calculus is a generalization of vector calculus, so for low-dimensional spaces, this translation is straightforward. For example, divergence is typically written as (⋆, d, ⋆). Scalar fields are typically interpreted as "0Forms", i.e. values assigned to vertices of a mesh.

We use the `@decapode` macro to interpret the equations. Here, we have equation 2 from Halfar:
```math
\frac{\partial h}{\partial t} = \frac{2}{n + 2} (\frac{\rho g}{A})^n \frac{\partial}{\partial x}(\frac{\partial h}{\partial x} |\frac{\partial h}{\partial x}| ^{n-1} h^{n+2}).
```
We'll change the term out front to Γ so we can demonstrate composition in a moment.

```
halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯(d(h)))^(n-1)) * avg₀₁(h^(n+2)))
end

to_graphviz(halfar_eq2)
```

And here, a formulation of Glen's law from J.W. Glen's 1958 "The flow law of ice".
```
glens_law = @decapode begin
  #Γ::Form0
  Γ::Form1
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

to_graphviz(glens_law)
```

# Composing models

We can use operadic composition to specify how our models come together. In this example, we have two Decapodes, and two quantities that are shared between them.
```
ice_dynamics_composition_diagram = @relation () begin
  dynamics(n,Γ)
  stress(Γ,n)
end

to_graphviz(ice_dynamics_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")
```

To a apply a composition, we specify which Decapodes to plug into those boxes, and what each calls the corresponding shared variables internally.

```
ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:h,:Γ,:n]),
  Open(glens_law, [:Γ,:n])])

ice_dynamics = apex(ice_dynamics_cospan)
to_graphviz(ice_dynamics)
```

# Provide a semantics

To interpret our composed Decapode, we need to specify what Discrete Exterior Calculus to interpret our quantities in. Let's choose the 1D Discrete Exterior Calculus:
```
ice_dynamics1D = expand_operators(ice_dynamics)
infer_types!(ice_dynamics1D, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(ice_dynamics1D, op1_res_rules_1D, op2_res_rules_1D)

to_graphviz(ice_dynamics1D)
```

# Define a mesh

We'll need a mesh to simulate on. Since this is a 1D mesh, we can go ahead and make one right now:
```
# This is a 1D mesh, consisting of edges and vertices.
s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
# 20 hundred vertices along a line, connected by edges.
add_vertices!(s′, 20, point=Point2D.(range(0, 10_000, length=20), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())
```

# Define input data

We need initial conditions to use for our simulation.
```
n = 3
ρ = 910
g = 9.8
A = 1e-16

# Ice height is a primal 0-form, with values at vertices.
# We choose a distribution that obeys the shallow height and shallow slope conditions.
h₀ = map(point(s′)) do (x,_)
  ((7072-((x-5000)^2))/9e3+2777)/2777e-1
end

# Visualize initial conditions for ice sheet height.
lines(map(x -> x[1], point(s′)), h₀, linewidth=5)
```

We need to tell our Decapode which data maps to which symbols. We can wrap up our data like so:
```
u₀ = construct(PhysicsState, [VectorForm(h₀)], Float64[], [:h])
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)
```

# Define functions

In order to solve our equations, we will need numerical linear operators that give meaning to our symbolic operators. In the DEC, there are a handful of operators that one uses to construct all the usual vector calculus operations, namely: ♯, ♭, ∧, d, ⋆. The CombinatorialSpaces.jl library specifies many of these for us.

```
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :♯ => x -> begin
      # This is an implementation of the "sharp" operator from the exterior
      # calculus, which takes co-vector fields to vector fields.
      # This could be up-streamed to the CombinatorialSpaces.jl library. (i.e.
      # this operation is not bespoke to this simulation.)
      e_vecs = map(edges(sd)) do e
        point(sd, sd[e, :∂v0]) - point(sd, sd[e, :∂v1])
      end
      neighbors = map(vertices(sd)) do v
        union(incident(sd, v, :∂v0), incident(sd, v, :∂v1))
      end
      n_vecs = map(neighbors) do es
        [e_vecs[e] for e in es]
      end
      map(neighbors, n_vecs) do es, nvs
        sum([nv*norm(nv)*x[e] for (e,nv) in zip(es,nvs)]) / sum(norm.(nvs))
      end
    end
    :mag => x -> begin
      norm.(x)
    end
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(s)
          append!(J, [s[e,:∂v0],s[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      avg_mat = sparse(I,J,V)
      avg_mat * x
    end
    :^ => (x,y) -> x .^ y
    :* => (x,y) -> x .* y
    :show => x -> begin
      @show x
      x
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

# Generate the simulation

Now, we have everything we need to generate our simulation:

```
sim = eval(gensim(ice_dynamics1D, dimension=1))
fₘ = sim(s, generate)
```

# Pre-compile and run

The first time that you run a function, Julia will pre-compile it, so that later runs will be fast. We'll solve our simulation for a short time span, to trigger this pre-compilation, and then run it.

```
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

tₑ = 8e3

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
```

We can save our solution file in case we want to examine its contents when this Julia session ends.

```
@save "ice_dynamics1D.jld2" soln
```

# Visualize

Let's examine the final conditions:

```
fig,ax,ob = lines(map(x -> x[1], point(s′)), findnode(soln(tₑ), :h), linewidth=5)
ylims!(ax, extrema(h₀))
display(fig)
```

We see that our distribution converges to a more uniform ice height across our domain, which matches our physical intuition.

Let's create a GIF to examine an animation of these dynamics:

```
# Create a gif
begin
  frames = 100
  fig, ax, ob = lines(map(x -> x[1], point(s′)), findnode(soln(0), :h))
  ylims!(ax, extrema(h₀))
  record(fig, "ice_dynamics1D.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    lines!(map(x -> x[1], point(s′)), findnode(soln(t), :h))
  end
end
```

![IceDynamics1D]("./ice_dynamics1D.gif")


# 2D Re-interpretation

The first, one-dimensional, semantics that we provided to our Decapode restricted the kinds of glacial sheets that we could model. (i.e. We could only look at glacial sheets which were constant along y). We can give our Decapode an alternate semantics, as some physics on a 2-dimensional manifold.

Note that for these physics, we make no adjustments to the underlying "dimension-agnostic" Decapode, we only provide a different set of rules for inferring what the type of each quantity is.

```
ice_dynamics2D = expand_operators(ice_dynamics)
infer_types!(ice_dynamics2D)
resolve_overloads!(ice_dynamics2D)

to_graphviz(ice_dynamics2D)
```

# Store as JSON

We quickly demonstrate how to serialize a Decapode to JSON and read it back in:
```
write_json_acset(ice_dynamics2D, "ice_dynamics2D.json")
```

```
# When reading back in, we specify that all attributes are "Symbol"s.
ice_dynamics2 = read_json_acset(SummationDecapode{Symbol,Symbol,Symbol}, "ice_dynamics2D.json")
# Or, you could choose to interpret the data as "String"s.
ice_dynamics3 = read_json_acset(SummationDecapode{String,String,String}, "ice_dynamics2D.json")

to_graphviz(ice_dynamics3)
```

# Define our mesh

```
include("../../grid_meshes.jl")
s′ = triangulated_grid(10_000,10_000,800,800,Point3D)
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())
wireframe(s)
```

# Define our input data

```
n = 3
ρ = 910
g = 9.8
A = 1e-16

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s′)) do (x,y)
  (7072-((x-5000)^2 + (y-5000)^2)^(1/2))/9e3+10
end

# Visualize initial condition for ice sheet height.
mesh(s′, color=h₀, colormap=:jet)
```

```
u₀ = construct(PhysicsState, [VectorForm(h₀)], Float64[], [:h])
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)
```

# Define our functions

```
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :♯ => x -> begin
      ♯(sd, EForm(x))
    end
    :mag => x -> begin
      norm.(x)
    end
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(s)
          append!(J, [s[e,:∂v0],s[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      avg_mat = sparse(I,J,V)
      avg_mat * x
    end
    :^ => (x,y) -> x .^ y
    :* => (x,y) -> x .* y
    :show => x -> begin
      @show x
      @show length(x)
      x
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

# Generate simulation

```
sim = eval(gensim(ice_dynamics2D, dimension=2))
fₘ = sim(s, generate)
```

# Pre-compile and run

```
@info("Precompiling Solver")
# We run for a short timespan to pre-compile.
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
```

```
tₑ = 5e13

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
```

```
@save "ice_dynamics2D.jld2" soln
```

# Visualize

```
# Final conditions:
mesh(s′, color=findnode(soln(tₑ), :h), colormap=:jet, colorrange=extrema(findnode(soln(0), :h)))
```

```
begin
  frames = 100
  fig, ax, ob = GLMakie.mesh(s′, color=findnode(soln(0), :h), colormap=:jet, colorrange=extrema(findnode(soln(0), :h)))
  Colorbar(fig[1,2], ob)
  record(fig, "ice_dynamics2D.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    ob.color = findnode(soln(t), :h)
  end
end
```

![IceDynamics2D]("./ice_dynamics2D.gif")

# 2-Manifold in 3D

We note that just because our physics is happening on a 2-manifold, (a surface), this doesn't restrict us to the 2D plane. In fact, we can "embed" our 2-manifold in 3D space to simulate a glacial sheets spread across the globe.

```
s′ = loadmesh(Icosphere(3, 10_000))
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())
wireframe(s)

```

```
n = 3
ρ = 910
g = 9.8
A = 1e-16

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s′)) do (x,y,z)
  (z*z)/(10_000*10_000)
end

# Visualize initial condition for ice sheet height.
# There is lots of ice at the poles, and no ice at the equator.
mesh(s′, color=h₀, colormap=:jet)
```

```
u₀ = construct(PhysicsState, [VectorForm(h₀)], Float64[], [:h])
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)
```

```
sim = eval(gensim(ice_dynamics2D, dimension=2))
fₘ = sim(s, generate)
```

```
# For brevity's sake, we'll skip the pre-compilation cell.

tₑ = 1e11

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

# Compare the extrema of the initial and final conditions of ice height.
extrema(findnode(soln(0), :h)), extrema(findnode(soln(tₑ), :h))
```

```
mesh(s′, color=findnode(soln(tₑ), :h), colormap=:jet, colorrange=extrema(findnode(soln(0), :h)))
```

```
begin
  frames = 200
  fig, ax, ob = GLMakie.mesh(s′, color=findnode(soln(0), :h), colormap=:jet, colorrange=extrema(findnode(soln(0), :h)))
  Colorbar(fig[1,2], ob)
  # These particular initial conditions diffuse quite quickly, so let's just look at
  # the first moments of those dynamics.
  record(fig, "ice_dynamics2D_sphere.gif", range(0.0, tₑ/64; length=frames); framerate = 20) do t
    ob.color = findnode(soln(t), :h)
  end
end
```

![IceDynamics2DSphere]("./ice_dynamics2D_sphere.gif")
