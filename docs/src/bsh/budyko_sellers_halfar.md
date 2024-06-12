# Budko-Sellers-Halfar

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..","docinfo.jl"))
info = DocInfo.Info()
```

In this example, we will compose the Budyko-Sellers 1D energy balance model of the Earth's surface temperature with the Halfar model of glacial dynamics. Note that each of these components models is itself a composition of smaller physical models. In this walkthrough, we will compose them together using the same techniques.

``` @example DEC
# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics: Point2
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays
Point2D = Point2{Float64};
nothing # hide
```

We have defined the [Halfar ice model](../ice_dynamics/ice_dynamics.md) in other docs pages, and so will quickly define it here.

``` @example DEC
halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ  * d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end

glens_law = @decapode begin
  Γ::Form1
  A::Form1
  (ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

ice_dynamics_composition_diagram = @relation () begin
  dynamics(Γ,n)
  stress(Γ,n)
end

ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:Γ,:n]),
   Open(glens_law, [:Γ,:n])])
halfar = apex(ice_dynamics_cospan)

to_graphviz(halfar, verbose=false)
```

We will introduce the Budyko-Sellers energy balance model in more detail. First, let's define the composite physics. We will visualize them all in a single diagram without any composition at first:

``` @example DEC
energy_balance = @decapode begin
  (Tₛ, ASR, OLR, HT)::Form0
  (C)::Constant

  Tₛ̇ == ∂ₜ(Tₛ) 

  Tₛ̇ == (ASR - OLR + HT) ./ C
end

absorbed_shortwave_radiation = @decapode begin
  (Q, ASR)::Form0
  α::Constant

  ASR == (1 .- α) .* Q
end

outgoing_longwave_radiation = @decapode begin
  (Tₛ, OLR)::Form0
  (A,B)::Constant

  OLR == A .+ (B .* Tₛ)
end

heat_transfer = @decapode begin
  (HT, Tₛ)::Form0
  (D,cosϕᵖ,cosϕᵈ)::Constant

  HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
end

insolation = @decapode begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end

to_graphviz(oplus([energy_balance, absorbed_shortwave_radiation, outgoing_longwave_radiation, heat_transfer, insolation]), directed=false)
```

Now let's compose the Budyko-Sellers model:

``` @example DEC
budyko_sellers_composition_diagram = @relation () begin
  energy(Tₛ, ASR, OLR, HT)
  absorbed_radiation(Q, ASR)
  outgoing_radiation(Tₛ, OLR)
  diffusion(Tₛ, HT, cosϕᵖ)
  insolation(Q, cosϕᵖ)
end

draw_composition(budyko_sellers_composition_diagram)
```

``` @example DEC
budyko_sellers_cospan = oapply(budyko_sellers_composition_diagram,
  [Open(energy_balance,               [:Tₛ, :ASR, :OLR, :HT]),
   Open(absorbed_shortwave_radiation, [:Q, :ASR]),
   Open(outgoing_longwave_radiation,  [:Tₛ, :OLR]),
   Open(heat_transfer,                [:Tₛ, :HT, :cosϕᵖ]),
   Open(insolation,                   [:Q, :cosϕᵖ])])

budyko_sellers = apex(budyko_sellers_cospan)

# Save this Decapode as a JSON file
write_json_acset(budyko_sellers, "budyko_sellers.json") 

to_graphviz(budyko_sellers, verbose=false)
```

## Warming

We need to specify physically what it means for these two terms to interact. We will say that ice will diffuse faster as temperature increases, and will pick some coefficients that demonstrate interesting dynamics on short timescales.

``` @example DEC
warming = @decapode begin
  Tₛ::Form0
  A::Form1

  A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.65e7)

end

to_graphviz(warming)
```

## Composition

Observe that Decapodes composition is hierarchical. This composition technique is the same as that used in composing each of the Budyko-Sellers and Halfar models.

``` @example DEC
budyko_sellers_halfar_composition_diagram = @relation () begin
  budyko_sellers(Tₛ)
  warming(A, Tₛ)
  halfar(A)
end

draw_composition(budyko_sellers_halfar_composition_diagram)
```

We apply a composition by plugging in a Decapode for each component. We also specify the internal name of the variables to be used in combining.

``` @example DEC
budyko_sellers_halfar_cospan = oapply(budyko_sellers_halfar_composition_diagram,
  [Open(budyko_sellers, [:Tₛ]),
   Open(warming,        [:A, :Tₛ]),
   Open(halfar,         [:stress_A])])
budyko_sellers_halfar = apex(budyko_sellers_halfar_cospan)

to_graphviz(budyko_sellers_halfar)
```

We can perform type inference to determine what kind of differential form each of our variables are. This is done automatically with the `dimension=1` keyword given to `gensim`, but we will do it in-place for demonstration purposes.

``` @example DEC
budyko_sellers_halfar = expand_operators(budyko_sellers_halfar)
infer_types!(budyko_sellers_halfar, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(budyko_sellers_halfar, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(budyko_sellers_halfar)
```

## Defining the mesh

These dynamics will occur on a 1-D manifold (a line). Points near +-π/2 will represent points near the North/ South poles. Points near 0 represent those at the equator.

``` @example DEC
s = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s, 100, point=Point2D.(range(-π/2 + π/32, π/2 - π/32, length=100), 0))
add_edges!(s, 1:nv(s)-1, 2:nv(s))
orient!(s)
sd = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s)
subdivide_duals!(sd, Circumcenter())
```

## Define input data

We need to supply initial conditions to our model. We will use synthetic data here.

``` @example DEC
# This is a primal 0-form, with values at vertices.
cosϕᵖ = map(x -> cos(x[1]), point(s))

# This is a dual 0-form, with values at edge centers.
cosϕᵈ = map(edges(s)) do e
  (cos(point(s, src(s, e))[1]) + cos(point(s, tgt(s, e))[1])) / 2
end

α₀ = 0.354
α₂ = 0.25
α = map(point(s)) do ϕ
  α₀ + α₂*((1/2)*(3*ϕ[1]^2 - 1))
end
A = 210
B = 2
f = 0.70
ρ = 1025
cw = 4186
H = 70
C = map(point(s)) do ϕ
  f * ρ * cw * H
end
D = 0.6

# Isothermal initial conditions:
Tₛ₀ = map(point(s)) do ϕ
  15
end

# Visualize initial condition for temperature.
lines(map(x -> x[1], point(s)), Tₛ₀)
```

``` @example DEC
n = 3
ρ = 910
g = 9.8

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s)) do (x,_)
  (((x)^2)+2.5) / 1e3
end

# Visualize initial condition for ice sheet height.
lines(map(x -> x[1], point(s)), h₀)
```

``` @example DEC
# Store these values to be passed to the solver.
u₀ = ComponentArray(Tₛ=Tₛ₀, halfar_dynamics_h=h₀)

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

## Symbols to functions

The symbols along edges in our Decapode must be mapped to executable functions. In the Discrete Exterior Calculus, all our operators are defined as relations between points, lines, and triangles on meshes known as simplicial sets. Thus, DEC operators are re-usable across any simplicial set.

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

## Simulation generation

From our Decapode, we automatically generate a finite difference method solver that performs explicit time-stepping to solve our system of multiphysics equations.

``` @example DEC
sim = eval(gensim(budyko_sellers_halfar, dimension=1))
fₘ = sim(sd, generate)
```

## Run simulation

We wrap our simulator and initial conditions and solve them with the stability-detection and time-stepping methods provided by DifferentialEquations.jl .

``` @example DEC
tₑ = 1e6

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
```

We can save the [solution file](budyko_sellers_halfar.jld2) to examine later.

``` @example DEC
@save "budyko_sellers_halfar.jld2" soln
```

## Visualize

Quickly examine the final conditions for temperature:

``` @example DEC
lines(map(x -> x[1], point(s)), soln(tₑ).Tₛ)
```

Quickly examine the final conditions for ice height:

``` @example DEC
lines(map(x -> x[1], point(s)), soln(tₑ).halfar_dynamics_h)
```

``` @setup DEC
begin
# Initial frame
frames = 100
fig = Figure()
ax1 = CairoMakie.Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s))))
ylims!(ax1, extrema(soln(tₑ).Tₛ))
ax1.xlabel = "Line plot of temperature from North to South pole, every $(tₑ/frames) time units"
Label(fig[1,1,Top()], "Surface temperature, Tₛ, [C°]")

# Animation
record(fig, "budyko_sellers_halfar_T.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s)), soln(t).Tₛ)
end
end

begin
# Initial frame
frames = 100
fig = Figure()
ax1 = CairoMakie.Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s))))
ylims!(ax1, extrema(soln(tₑ).halfar_dynamics_h))
ax1.xlabel = "Line plot of temperature from North to South pole, every $(tₑ/frames) time units"
Label(fig[1,1,Top()], "Ice height, h")

# Animation
record(fig, "budyko_sellers_halfar_h.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s)), soln(t).halfar_dynamics_h)
end
end
```

![BSH_Temperature](budyko_sellers_halfar_T.gif)

![BSH_IceHeight](budyko_sellers_halfar_h.gif)

```@example INFO
DocInfo.get_report(info) # hide
```
