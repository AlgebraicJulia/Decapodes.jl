# Brusselator

```@setup INFO
include(joinpath(Base.@__DIR__, ".." , "..", "docinfo.jl"))
info = DocInfo.Info()
```
## Dependencies
```@example DEC
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis

using CombinatorialSpaces
using ComponentArrays
using DiagrammaticEquations
using Decapodes
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}
```

## The Model
```@example DEC
Brusselator = @decapode begin
  (U, V)::Form0
  U2V::Form0
  (U̇, V̇)::Form0

  (α)::Constant
  F::Parameter

  U2V == (U .* U) .* V

  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(V))
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end
```

## The Mesh
```@example DEC
s = triangulated_grid(1,1,0.008,0.008,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());
```

## Initial data
```@example DEC
U = map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end

V = map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end

F₁ = map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end

F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₁ : F₂)
```

## Generate the Simulation
```@example DEC
sim = evalsim(Brusselator)
fₘ = sim(sd, nothing, DiagonalHodge())

u₀ = ComponentArray(U=U, V=V)

tₑ = 11.5

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
```

## Visualize
```@example DEC

fig = Figure();
ax = Axis(fig[1,1])
mesh!(ax, s, color=soln(tₑ).U, colormap=:jet)

function save_dynamics(save_file_name)
  time = Observable(0.0)
  u = @lift(soln($time).U)
  f = Figure()
  ax = CairoMakie.Axis(f[1,1], title = @lift("Brusselator U Concentration at Time $($time)"))
  gmsh = mesh!(ax, s, color=u, colormap=:jet,
               colorrange=extrema(soln(tₑ).U))
  Colorbar(f[1,2], gmsh)
  timestamps = range(0, tₑ, step=1e-1)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("brusselator.gif")
```

![Brusselator_results_flat](brusselator.gif)

```@example INFO
DocInfo.get_report(info) # hide
```

