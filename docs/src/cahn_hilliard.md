# The Cahn-Hilliard Equations

For this example, Decapodes will model the Cahn-Hilliard equations. These equations describe the evolution of a binary fluid as its two phases seperate out into distinct domains.

## Formulating the Equations

We first load in our dependencies.

```@example DEC
# AlgebraicJulia Dependencies
using Catlab
using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces

# External Dependencies
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie

using GeometryBasics
Point3D = Point3{Float64};
nothing
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

In these equations, C will represent the concentration of the binary fluid, range from -1 to 1 to differentiate between different phase. Then we have D as a diffusion constant and γ as the length of the transition regions. This formulation of the Cahn-Hilliard equations was drawn from the Wikipedia page on the topic, found here: https://en.wikipedia.org/wiki/Cahn%E2%80%93Hilliard_equation.

## Loading the Data

We now generate the mesh information. We first run the equations on a triangulated grid.

```@example DEC
s = triangulated_grid(100, 100, 0.5, 0.5, Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter());

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect=1)
wf = wireframe!(ax, s)
save("CahnHillard_Rect.png", fig)
```

The Cahn-Hilliard equations start with random concentration values between -1 and 1. For both D and γ constants we choose 0.5.

```@example DEC
C = rand(Float64, nv(sd)) * 2 .- 1
u₀ = ComponentArray(C=C)

constants = (D = 0.5, γ = 0.5);
```

We'll now create the simulation code representing the Cahn-Hilliard equations. We pass nothing in the second arguement to sim since we have no custom functions to pass in.

```@example DEC
sim = eval(gensim(CahnHilliard))
fₘ = sim(sd, nothing, DiagonalHodge());
```

## Getting the Solution

Now that everything is set up and ready, we can solve the equations. We run the simulation for 240 time units, to see the long-term evolution of the fluid. Note we only save the solution at intervals of 0.1 time units in order to save memory.

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

create_gif(soln, "CahnHillard_Rect.gif")
```