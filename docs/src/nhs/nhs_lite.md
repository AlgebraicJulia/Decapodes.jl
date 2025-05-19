# Implement Oceananigans.jl's NonhydrostaticModel in the Discrete Exterior Calculus

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..", "docinfo.jl"))
info = DocInfo.Info()
```

Let's use Decapodes to implement the [NonhydrostaticModel](https://clima.github.io/OceananigansDocumentation/stable/physics/nonhydrostatic_model/) from Oceananigans.jl. We will take the opportunity to demonstrate how we can use our "algebra of model compositions" to encode certain guarantees on the models we generate. We will use the [2D Turbulence](https://clima.github.io/OceananigansDocumentation/stable/literated/two_dimensional_turbulence/) as a guiding example, and use only equations found in the Oceananigans docs to construct our model.

The full code that generated these results is available in a [julia script](nhs.jl).

```@example DEC
# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations

# External Dependencies
using CairoMakie
using ComponentArrays
using Downloads
using GeometryBasics: Point3
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
Point3D = Point3{Float64};
nothing # hide
```

## Specify our models

This is [Equation 1: "The momentum conservation equation"](https://clima.github.io/OceananigansDocumentation/stable/physics/nonhydrostatic_model/#The-momentum-conservation-equation). This is the first formulation of mutual advection (of `v` along `V`, and `V` along `v`) that we could find in the exterior calculus.

```@example DEC
momentum =  @decapode begin
  (v,V)::DualForm1
  f::Form0
  uˢ::DualForm1
  ∂tuˢ::DualForm1
  p::DualForm0
  b::DualForm0
  ĝ::DualForm1
  Fᵥ::DualForm1
  StressDivergence::DualForm1

  ∂ₜ(v) ==
    -ℒ₁(v,v) + 0.5*d(ι₁₁(v,v)) -
     d(ι₁₁(v,V)) + ι₁₂(v,d(V)) + ι₁₂(V,d(v)) -
     (f - ∘(d,⋆)(uˢ)) ∧ᵖᵈ₀₁ v -
     d(p) +
     b ∧ᵈᵈ₀₁ ĝ -
     StressDivergence +
     ∂tuˢ +
     Fᵥ
end
to_graphviz(momentum)
```

Why did we write "StressDivergence" instead of ∇⋅τ, as in the linked equation? According to [this docs page](https://clima.github.io/OceananigansDocumentation/stable/physics/turbulence_closures/), the user makes a selection of what model to insert in place of the term ∇⋅τ. For example, in [the isotropic case](https://clima.github.io/OceananigansDocumentation/stable/physics/turbulence_closures/#Constant-isotropic-diffusivity), Oceananigans.jl replaces this term with: ∇⋅τ = *ν*Δv. Thus, we write StressDivergence, and replace this term with a choice of "turbulence closure" model. Using the "constant isotropic diffusivity" case, we can operate purely in terms of scalar-valued forms.

This is [Equation 2: "The tracer conservation equation"](https://clima.github.io/OceananigansDocumentation/stable/physics/nonhydrostatic_model/#The-tracer-conservation-equation).

```@example DEC
tracer_conservation = @decapode begin
  (c,C,F,FluxDivergence)::DualForm0
  (v,V)::DualForm1

  ∂ₜ(c) ==
    -1*ι₁₁(v,d(c)) -
    ι₁₁(V,d(c)) -
    ι₁₁(v,d(C)) -
    FluxDivergence +
    F
end
to_graphviz(tracer_conservation)
```

This is [Equation 2: "Linear equation of state"](https://clima.github.io/OceananigansDocumentation/stable/physics/buoyancy_and_equations_of_state/#Linear-equation-of-state) of seawater buoyancy.

```@example DEC
equation_of_state = @decapode begin
  (b,T,S)::DualForm0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end
to_graphviz(equation_of_state)
```

This is [Equation 2: "Constant isotropic diffusivity"](https://clima.github.io/OceananigansDocumentation/stable/physics/turbulence_closures/#Constant-isotropic-diffusivity).

```@example DEC
isotropic_diffusivity = @decapode begin
  v::DualForm1
  c::DualForm0
  StressDivergence::DualForm1
  FluxDivergence::DualForm0
  (κ,nu)::Constant

  StressDivergence == nu*Δᵈ₁(v)
  FluxDivergence == κ*Δᵈ₀(c)
end
to_graphviz(isotropic_diffusivity)
```

## Compatibility Guarantees via Operadic Composition

Decapodes composition is formally known as an "operad algebra". That means that we don't have to encode our composition in a single undirected wiring diagram (UWD) and then apply it. Rather, we can define several UWDs, compose those, and then apply those. Of course, since the output of oapply is another Decapode, we could perform an intermediate oapply, if that is convenient.

Besides it being convenient to break apart large UWDs into component UWDs, this hierarchical composition can enforce rules on our physical quantities.

For example:

1. We want all the tracers (salinity, temperature, etc.) in our physics to obey the same conservation equation.
2. We want them to obey the same "turbulence closure", which affects their flux-divergence term.
3. At the same time, a choice of turbulence closure doesn't just affect (each of) the flux-divergence terms, it also constrains which stress-divergence is physically valid in the momentum equation.

We will use our operad algebra to guarantee model compatibility and physical consistency, guarantees that would be burdensome to fit into a one-off type system.

Here, we specify the equations that any tracer obeys:

```@example DEC
tracer_composition = @relation () begin
  # "The turbulence closure selected by the user determines the form of ... diffusive flux divergence"
  turbulence(FD,v,c)

  continuity(FD,v,c)
end
draw_composition(tracer_composition)
```

Let's "lock in" isotropic diffusivity by doing an intermediate oapply.

```@example DEC
isotropic_tracer = apex(oapply(tracer_composition, [
  Open(isotropic_diffusivity, [:FluxDivergence, :v, :c]),
  Open(tracer_conservation,   [:FluxDivergence, :v, :c])]))
to_graphviz(isotropic_tracer)
```

Let's use this building-block tracer physics at the next level. The quotes that appear in this composition diagram appear directly in the Oceananigans.jl docs.

```@example DEC
nonhydrostatic_composition = @relation () begin
  # "The turbulence closure selected by the user determines the form of stress divergence"
  #   => Note that the StressDivergence term, SD, is shared by momentum and all the tracers.
  momentum(V, v, b, SD)

  # "Both T and S obey the tracer conservation equation"
  #   => Temperature and Salinity both receive a copy of the tracer physics.
  temperature(V, v, T, SD, nu)
  salinity(V, v, S, SD, nu)

  # "Buoyancy is determined from a linear equation of state"
  #   => The b term in momentum is that described by the equation of state here.
  eos(b, T, S)
end
draw_composition(nonhydrostatic_composition)
```

```@example DEC
isotropic_nonhydrostatic_buoyancy = apex(oapply(nonhydrostatic_composition, [
  Open(momentum,          [:V, :v, :b, :StressDivergence]),
  Open(isotropic_tracer,  [:continuity_V, :v, :c, :turbulence_StressDivergence, :turbulence_nu]),
  Open(isotropic_tracer,  [:continuity_V, :v, :c, :turbulence_StressDivergence, :turbulence_nu]),
  Open(equation_of_state, [:b, :T, :S])]));
to_graphviz(isotropic_nonhydrostatic_buoyancy)
```

## Our Mesh

We execute these dynamics on the torus explicitly, instead of using a square with periodic boundary conditions.

```@example DEC
# This is a torus with resolution of its dual mesh similar to that
# used by Oceananigans (explicitly represented as a torus, not as a
# square with periodic boundary conditions!)
Downloads.download("https://cise.ufl.edu/~luke.morris/torus.obj", "torus.obj")
s = EmbeddedDeltaSet2D(joinpath(@__DIR__, "torus.obj"))
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Barycenter())
fig = Figure() 
ax = CairoMakie.Axis(fig[1,1], aspect=1) 
wf = wireframe!(ax, s; linewidth=1) 
save("NHS_mesh.png", fig) 
nothing # hide
```

!["NHS_torus"](NHS_mesh.png)
  
## Results

In the DEC, vorticity is encoded with `d⋆`, and speed can be encoded with `norm ♯`. We can use our operators from CombinatorialSpaces.jl to create our GIFs.

![Vorticity](vorticity.gif)

![Speed](speed.gif)

```@example INFO
DocInfo.get_report(info) # hide
```
