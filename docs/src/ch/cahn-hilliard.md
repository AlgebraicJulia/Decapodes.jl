# The Cahn-Hilliard Equation

```@setup DEC
using Dates
start_time = now()
```

For this example Decapodes will model the Cahn-Hilliard equation. This equation describes the evolution of a binary fluid as its two phases separate out into distinct domains.

!["Cahn Hilliard sample"](CahnHilliard_Final.jpg)

## Formulating the Equation

We first load in our dependencies.

```@example DEC
# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using Random
Point3D = Point3{Float64};
nothing #hide
```

We then proceed to describe our physics using Decapodes.

```@example DEC
CahnHilliard = @decapode begin
    C::Form0
    (D, γ)::Constant
    ∂ₜ(C) == D * Δ(C.^3 - C - γ * Δ(C))
end

to_graphviz(CahnHilliard)
```

In this equation C will represent the concentration of the binary fluid, ranging from -1 to 1 to differentiate between different phases. We also have a diffusion constant D and a constant γ whose square root is the length of the transition regions. This formulation of the Cahn-Hilliard equation was drawn from the Wikipedia page on the topic found [here](https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation).

## Loading the Data

We now generate the mesh information. We'll run the equation on a triangulated grid.

```@example DEC
s = triangulated_grid(100, 100, 0.8, 0.8, Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter());
fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect=1)
wf = wireframe!(ax, s; linewidth=1)
save("CahnHilliard_Rect.png", fig)
```

!["CahnHilliardRect"](CahnHilliard_Rect.png)

The Cahn-Hilliard equation starts with a random concentration holding values between -1 and 1. For both D and γ constants we choose 0.5.

```@example DEC
Random.seed!(0)

C = rand(Float64, nv(sd)) * 2 .- 1
u₀ = ComponentArray(C=C)
constants = (D = 0.5, γ = 0.5);
```

We'll now create the simulation code representing the Cahn-Hilliard equation. We pass nothing in the second argument to sim since we have no custom functions to pass in.

```@example DEC
sim = eval(gensim(CahnHilliard))
fₘ = sim(sd, nothing, DiagonalHodge());
```

## Getting the Solution

Now that everything is set up and ready, we can solve the equation. We run the simulation for 240 time units to see the long-term evolution of the fluid. Note we only save the solution at intervals of 0.1 time units in order to reduce the memory-footprint of the solve.

```@example DEC
tₑ = 240
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants)
soln = solve(prob, Tsit5(), saveat=0.1);
soln.retcode
```

And we can see the result as a gif.

```@example DEC
function create_gif(solution, file_name)
  frames = 200
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s, color=solution(0).C, colormap=:jet, colorrange=extrema(solution(0).C))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, file_name, range(0.0, tₑ; length=frames); framerate = 15) do t
    msh.color = solution(t).C
  end
end
create_gif(soln, "CahnHilliard_Rect.gif")
```

!["CahnHilliardRes"](CahnHilliard_Rect.gif)

```@example DEC
end_time = now() # hide
elapsed = end_time - start_time # hide
elapsed_sec = round(elapsed, Dates.Second(1)) # hide
@info "Documentation completed in $(elapsed_sec)." # hide
@info "Documentation built at $(end_time)." # hide
```
