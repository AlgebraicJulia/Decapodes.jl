# Halfar's model of glacial flow

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..", "docinfo.jl"))
info = DocInfo.Info()
```

Let's model glacial flow using a model of how ice height of a glacial sheet changes over time, from P. Halfar's 1981 paper: ["On the dynamics of the ice sheets"](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JC086iC11p11065).

``` @example DEC
# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics: Point2, Point3
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using Statistics
Point2D = Point2{Float64};
Point3D = Point3{Float64};
```

## Defining the models

The first step is to find a suitable equation for our model, and translate it into the Discrete Exterior Calculus. The Exterior Calculus is a generalization of vector calculus, so for low-dimensional spaces, this translation is straightforward. For example, divergence is typically written as (⋆, d, ⋆). Scalar fields are typically interpreted as "0Forms", i.e. values assigned to vertices of a mesh.

We use the `@decapode` macro to interpret the equations. Here, we have equation 2 from Halfar:

```math
\frac{\partial h}{\partial t} = \frac{2}{n + 2} (\frac{\rho g}{A})^n \frac{\partial}{\partial x}(\frac{\partial h}{\partial x} |\frac{\partial h}{\partial x}| ^{n-1} h^{n+2}).
```

We'll change the term out front to Γ so we can demonstrate composition in a moment.

In the exterior calculus, we could write the above equations like so:

```math
\partial_t(h) = \Gamma\quad \circ(\star, d, \star)(d(h)\quad \wedge \quad|d(h)^\sharp|^{n-1} \quad \wedge \quad (h^{n+2})).
```

``` @example DEC
halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form0
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == Γ * ∘(⋆, d, ⋆)(d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end

to_graphviz(halfar_eq2)
```

And here, a formulation of Glen's law from J.W. Glen's 1958 ["The flow law of ice"](http://go.owu.edu/~chjackso/Climate/papers/Glen_1958_The%20flow%20law%20of%20ice.pdf).

``` @example DEC
glens_law = @decapode begin
  Γ::Form0
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

to_graphviz(glens_law)
```

## Composing models

We can use [operadic composition](https://arxiv.org/abs/2105.12282) to specify how our models come together. In this example, we have two Decapodes, and two quantities that are shared between them.

``` @example DEC
ice_dynamics_composition_diagram = @relation () begin
  dynamics(Γ,n)
  stress(Γ,n)
end

draw_composition(ice_dynamics_composition_diagram)
```

To a apply a composition, we specify which Decapodes to plug into those boxes, and what each calls the corresponding shared variables internally.

``` @example DEC
ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:Γ,:n]),
  Open(glens_law, [:Γ,:n])])

ice_dynamics = apex(ice_dynamics_cospan)
to_graphviz(ice_dynamics)
```

## Provide a semantics

To interpret our composed Decapode, we need to specify what Discrete Exterior Calculus to interpret our quantities in. Let's choose the 1D Discrete Exterior Calculus:

``` @example DEC
ice_dynamics1D = expand_operators(ice_dynamics)
infer_types!(ice_dynamics1D, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(ice_dynamics1D, op1_res_rules_1D, op2_res_rules_1D)

to_graphviz(ice_dynamics1D)
```

## Define a mesh

We'll need a mesh to simulate on. Since this is a 1D mesh, we can go ahead and make one right now:

``` @example DEC
# This is an empty 1D mesh.
s = EmbeddedDeltaSet1D{Bool, Point2D}()

# 20 vertices along a line, connected by edges.
add_vertices!(s, 20, point=Point2D.(range(0, 10_000, length=20), 0))
add_edges!(s, 1:nv(s)-1, 2:nv(s))
orient!(s)

# The dual 1D mesh
sd = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s)
subdivide_duals!(sd, Circumcenter())
```

## Define input data

We need initial conditions to use for our simulation.

``` @example DEC
n = 3
ρ = 910
g = 9.8
A = 1e-16

# Ice height is a primal 0-form, with values at vertices.
# We choose a distribution that obeys the shallow height and shallow slope conditions.
h₀ = map(point(s)) do (x,_)
  ((7072-((x-5000)^2))/9e3+2777)/2777e-1
end

# Visualize initial conditions for ice sheet height.
lines(map(x -> x[1], point(s)), h₀, linewidth=5)
```

We need to tell our Decapode which data maps to which symbols. We can wrap up our data like so:

``` @example DEC
u₀ = ComponentArray(dynamics_h=h₀)

constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)
```

## Define functions

In order to solve our equations, we will need numerical linear operators that give meaning to our symbolic operators. In the DEC, there are a handful of operators that one uses to construct all the usual vector calculus operations, namely: ♯, ♭, ∧, d, ⋆. The CombinatorialSpaces.jl library specifies many of these for us.

``` @example DEC
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
    :mag => x -> norm.(x)
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

## Generate the simulation

Now, we have everything we need to generate our simulation:

``` @example DEC
sim = eval(gensim(ice_dynamics1D, dimension=1))
fₘ = sim(sd, generate)
```

## Pre-compile and run

The first time that you run a function, Julia will pre-compile it, so that later runs will be fast. We'll solve our simulation for a short time span, to trigger this pre-compilation, and then run it.

``` @example DEC
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

tₑ = 8_000

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
```

We can save our [solution file](ice_dynamics1D.jld2) in case we want to examine its contents when this Julia session ends.

``` @example DEC
@save "ice_dynamics1D.jld2" soln
```

## Visualize

Let's examine the final conditions:

``` @example DEC
fig,ax,ob = lines(map(x -> x[1], point(s)), soln(tₑ).dynamics_h, linewidth=5)
ylims!(ax, extrema(h₀))
fig
```

We see that our distribution converges to a more uniform ice height across our domain, which matches our physical intuition.

Let's create a GIF to examine an animation of these dynamics:

``` @example DEC
# Create a gif
begin
  frames = 100
  fig, ax, ob = lines(map(x -> x[1], point(s)), soln(0).dynamics_h)
  ylims!(ax, extrema(h₀))
  record(fig, "ice_dynamics1D.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    lines!(map(x -> x[1], point(s)), soln(t).dynamics_h)
  end
end
```

![IceDynamics1D](ice_dynamics1D.gif)

## 2D Re-interpretation

The first, one-dimensional, semantics that we provided to our Decapode restricted the kinds of glacial sheets that we could model. (i.e. We could only look at glacial sheets which were constant along y). We can give our Decapode an alternate semantics, as some physics on a 2-dimensional manifold.

Note that for these physics, we make no adjustments to the underlying "dimension-agnostic" Decapode, we only provide a different set of rules for inferring what the type of each quantity is.

``` @example DEC
ice_dynamics2D = expand_operators(ice_dynamics)
infer_types!(ice_dynamics2D)
resolve_overloads!(ice_dynamics2D)

to_graphviz(ice_dynamics2D)
```

## Store as JSON

We quickly demonstrate how to serialize a Decapode to JSON and read it back in:

``` @example DEC
write_json_acset(ice_dynamics2D, "ice_dynamics2D.json")
```

You can view the JSON file [here](ice_dynamics2D.json).

``` @example DEC
# When reading back in, we specify that all attributes are "Symbol"s.
ice_dynamics2 = read_json_acset(SummationDecapode{Symbol,Symbol,Symbol}, "ice_dynamics2D.json")
# Or, you could choose to interpret the data as "String"s.
ice_dynamics3 = read_json_acset(SummationDecapode{String,String,String}, "ice_dynamics2D.json")

to_graphviz(ice_dynamics3)
```

## Define our mesh

``` @example DEC
s = triangulated_grid(10_000,10_000,800,800,Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Barycenter())

fig = Figure()
ax = CairoMakie.Axis(fig[1,1])
wf = wireframe!(ax, s)
fig
```

## Define our input data

``` @example DEC
n = 3
ρ = 910
g = 9.8
A = 1e-16

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s)) do (x,y)
  (7072-((x-5000)^2 + (y-5000)^2)^(1/2))/9e3+10
end

# Visualize initial condition for ice sheet height.
mesh(s, color=h₀, colormap=:jet)
fig = Figure()
ax = CairoMakie.Axis(fig[1,1])
msh = mesh!(ax, s, color=h₀, colormap=:jet)
Colorbar(fig[1,2], msh)
fig
```

``` @example DEC
u₀ = ComponentArray(dynamics_h=h₀)

constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)
```

## Define our functions

``` @example DEC
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :♯ => begin
      sharp_mat = ♯_mat(sd, AltPPSharp())
      x -> sharp_mat * x
    end
    :mag => x -> norm.(x)
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

## Generate simulation

``` @example DEC
sim = eval(gensim(ice_dynamics2D, dimension=2))
fₘ = sim(sd, generate)
```

## Pre-compile and run 2D

``` @example DEC
@info("Precompiling Solver")
# We run for a short timespan to pre-compile.
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
```

``` @example DEC
tₑ = 5e13

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
```

``` @example DEC
@save "ice_dynamics2D.jld2" soln
```

## Visualize 2D

``` @example DEC
# Final conditions:
fig = Figure()
ax = CairoMakie.Axis(fig[1,1])
msh = mesh!(ax, s, color=soln(tₑ).dynamics_h, colormap=:jet, colorrange=extrema(soln(0).dynamics_h))
Colorbar(fig[1,2], msh)
fig
```

``` @example DEC
begin
  frames = 100
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s, color=soln(0).dynamics_h, colormap=:jet, colorrange=extrema(soln(0).dynamics_h))
  Colorbar(fig[1,2], msh)
  record(fig, "ice_dynamics2D.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    msh.color = soln(t).dynamics_h
  end
end
```

![IceDynamics2D](ice_dynamics2D.gif)

## 2-Manifold in 3D

We note that just because our physics is happening on a 2-manifold, (a surface), this doesn't restrict us to the 2D plane. In fact, we can "embed" our 2-manifold in 3D space to simulate a glacial sheets spread across the globe.

``` @example DEC
s = loadmesh(Icosphere(3, 10_000))
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Barycenter())
wireframe(sd)
```

``` @example DEC
n = 3
ρ = 910
g = 9.8
A = 1e-16

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s)) do (x,y,z)
  (z*z)/(10_000*10_000)
end

# Visualize initial condition for ice sheet height.
# There is lots of ice at the poles, and no ice at the equator.
fig = Figure()
ax = LScene(fig[1,1], scenekw=(lights=[],))
msh = CairoMakie.mesh!(ax, s, color=h₀, colormap=:jet)
Colorbar(fig[1,2], msh)
fig
```

``` @example DEC
u₀ = ComponentArray(dynamics_h=h₀)

constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)
```

``` @example DEC
sim = eval(gensim(ice_dynamics2D, dimension=2))
fₘ = sim(sd, generate)
```

For brevity's sake, we'll skip the pre-compilation cell.

``` @example DEC
tₑ = 5e25

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

# Compare the extrema of the initial and final conditions of ice height.
extrema(soln(0).dynamics_h), extrema(soln(tₑ).dynamics_h)
```

``` @example DEC
fig = Figure()
ax = LScene(fig[1,1], scenekw=(lights=[],))
msh = CairoMakie.mesh!(ax, s, color=soln(tₑ).dynamics_h, colormap=:jet, colorrange=extrema(soln(0).dynamics_h))
Colorbar(fig[1,2], msh)
fig
```

``` @example DEC
begin
  frames = 200
  fig = Figure()
  ax = LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.mesh!(ax, s, color=soln(0).dynamics_h, colormap=:jet, colorrange=extrema(soln(0).dynamics_h))

  Colorbar(fig[1,2], msh)
  # These particular initial conditions diffuse quite quickly, so let's just look at
  # the first moments of those dynamics.
  record(fig, "ice_dynamics2D_sphere.gif", range(0.0, tₑ/64; length=frames); framerate = 20) do t
    msh.color = soln(t).dynamics_h
  end
end
```

![IceDynamics2DSphere](ice_dynamics2D_sphere.gif)

```@example INFO
DocInfo.get_report(info) # hide
```
