# Budko-Sellers-Halfar

In this example, we will compose the Budyko-Sellers 1D energy balance model of the Earth's surface temperature with the Halfar model of glacial dynamics. Note that each of these components models is itself a composition of smaller physical models. In this walkthrough, we will compose them together using the same techniques.

``` @example DEC
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
using CairoMakie
using GeometryBasics: Point2
Point2D = Point2{Float64};
```

We defined the Budyko-Sellers and Halfar models in example scripts (soon to be turned into Docs pages) in the `examples/climate` folder of the main repository. We recall them here.

``` @setup DEC
include("../../examples/climate/budyko_sellers.jl")
```

``` @setup DEC
include("../../examples/climate/shallow_ice.jl")
```

``` @example DEC
budyko_sellers = apex(budyko_sellers_cospan)
halfar = apex(ice_dynamics_cospan)
true # hide
```

## Budyko-Sellers

``` @example DEC
to_graphviz(budyko_sellers)
```

## Halfar

``` @example DEC
to_graphviz(halfar)
```

## Warming

This is a formula that computes `A` for use in the Halfar glacial dynamics, given `T` from the Budyko-Sellers model.

``` @example DEC
# Tₛ(ϕ,t) := Surface temperature
# A(ϕ) := Longwave emissions at 0°C
warming = @decapode begin
  (Tₛ)::Form0
  (A)::Form1

  A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.65e7)

end
to_graphviz(warming)
```

## Composition

Observe that this composition technique is the same as that used in composing each of the Budyko-Sellers and Halfar models.

``` @example DEC
budyko_sellers_halfar_composition_diagram = @relation () begin
  budyko_sellers(Tₛ)

  warming(A, Tₛ)

  halfar(A)
end
to_graphviz(budyko_sellers_halfar_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")
```

We apply a composition by plugging in a Decapode for each component. We also specify the internal name of the variables to be used in combining.

``` @example DEC
budyko_sellers_halfar_cospan = oapply(budyko_sellers_halfar_composition_diagram,
  [Open(budyko_sellers, [:Tₛ]),
   Open(warming, [:A, :Tₛ]),
   Open(halfar, [:stress_A])])
budyko_sellers_halfar = apex(budyko_sellers_halfar_cospan)
to_graphviz(budyko_sellers_halfar)
```

We can perform type inference to determine what kind of differential form each of our variables are.

``` @example DEC
budyko_sellers_halfar = expand_operators(budyko_sellers_halfar)
infer_types!(budyko_sellers_halfar, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(budyko_sellers_halfar, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(budyko_sellers_halfar)
```

## Defining the mesh

These dynamics will occur on a 1-D manifold (a line). Points near +-π/2 will represent points near the North/ South poles. Points near 0 represent those at the equator.

``` @example DEC
s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
#add_vertices!(s′, 30, point=Point2D.(range(-π/2 + π/32, π/2 - π/32, length=30), 0))
add_vertices!(s′, 100, point=Point2D.(range(-π/2 + π/32, π/2 - π/32, length=100), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())
```

## Define input data

We need to supply initial conditions to our model. We create synthetic data here, although one may imagine that they could source this from their data repo of choice.

``` @example DEC
# This is a primal 0-form, with values at vertices.
cosϕᵖ = map(x -> cos(x[1]), point(s′))
# This is a dual 0-form, with values at edge centers.
cosϕᵈ = map(edges(s′)) do e
  (cos(point(s′, src(s′, e))[1]) + cos(point(s′, tgt(s′, e))[1])) / 2
end

α₀ = 0.354
α₂ = 0.25
α = map(point(s′)) do ϕ
  α₀ + α₂*((1/2)*(3*ϕ[1]^2 - 1))
end
A = 210
B = 2
f = 0.70
ρ = 1025
cw = 4186
H = 70
C = map(point(s′)) do ϕ
  f * ρ * cw * H
end
D = 0.6

# Isothermal initial conditions:
Tₛ₀ = map(point(s′)) do ϕ
  15
end

# Visualize initial condition for temperature.
lines(map(x -> x[1], point(s′)), Tₛ₀)
```

``` @example DEC
n = 3
ρ = 910
g = 9.8

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s′)) do (x,_)
  (((x)^2)+2.5) / 1e3
end

# Visualize initial condition for ice sheet height.
lines(map(x -> x[1], point(s′)), h₀)
```

``` @example DEC
# Store these values to be passed to the solver.
u₀ = construct(PhysicsState, [
  VectorForm(Tₛ₀),
  VectorForm(h₀),
  ], Float64[], [
  :Tₛ,
  :halfar_h])
constants_and_parameters = (
  budyko_sellers_absorbed_radiation_α = α,
  budyko_sellers_outgoing_radiation_A = A,
  budyko_sellers_outgoing_radiation_B = B,
  budyko_sellers_energy_C = C,
  budyko_sellers_diffusion_D = D,
  budyko_sellers_cosϕᵖ = cosϕᵖ,
  budyko_sellers_diffusion_cosϕᵈ = cosϕᵈ,
  halfar_n = n,
  halfar_stress_ρ = ρ,
  halfar_stress_g = g)
```

# Symbols to functions

The symbols along edges in our Decapode must be mapped to executable functions. In the Discrete Exterior Calculus, all our operators are defined as relations bewteen points, lines, and triangles on meshes known as simplicial sets. Thus, DEC operators are re-usable across any simplicial set.

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

## Simulation generation

From our Decapode, we automatically generate a finite difference method solver that performs explicit time-stepping to solve our system of multiphysics equations.

``` @example DEC
sim = eval(gensim(budyko_sellers_halfar, dimension=1))
fₘ = sim(s, generate)
```

## Run simulation

We wrap our simulator and initial conditions and solve them with the stability-detection and time-stepping methods provided by DifferentialEquations.jl .

``` @example DEC
tₑ = 1e6

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
```

We can save the solution file to examine later.

``` @example DEC
@save "budyko_sellers_halfar.jld2" soln
```

## Visualize

Quickly examine the final conditions for temperature.
``` @example DEC
lines(map(x -> x[1], point(s′)), findnode(soln(tₑ), :Tₛ))
```

Quickly examine the final conditions for ice height.
``` @example DEC
lines(map(x -> x[1], point(s′)), findnode(soln(tₑ), :halfar_h))
```

Create animated GIFs of the temperature and ice height dynamics.
``` @example DEC
begin
# Initial frame
frames = 100
fig = Figure(resolution = (800, 800))
ax1 = Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(findnode(soln(tₑ), :Tₛ)))
Label(fig[1,1,Top()], "Surface temperature, Tₛ, [C°]")
Label(fig[2,1,Top()], "Line plot of temperature from North to South pole, every $(tₑ/frames) time units")

# Animation
record(fig, "budyko_sellers_halfar_T.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), findnode(soln(t), :Tₛ))
end
end

begin
# Initial frame
frames = 100
fig = Figure(resolution = (800, 800))
ax1 = Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(findnode(soln(tₑ), :halfar_h)))
Label(fig[1,1,Top()], "Ice height, h")
Label(fig[2,1,Top()], "Line plot of ice height from North to South pole, every $(tₑ/frames) time units")

# Animation
record(fig, "budyko_sellers_halfar_h.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), findnode(soln(t), :halfar_h))
end
end
```

![BSH_Temperature](budyko_sellers_halfar_T.gif)

![BSH_IceHeight](budyko_sellers_halfar_h.gif)
