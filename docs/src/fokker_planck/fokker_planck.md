# Fokker-Planck

``` @example DEC
using Catlab, CombinatorialSpaces, Decapodes, DiagrammaticEquations
using CairoMakie, ComponentArrays, LinearAlgebra, MLStyle, ComponentArrays
using OrdinaryDiffEq
using GeometryBasics: Point3
Point3D = Point3{Float64}
using Arpack
```

Let's specify physics
``` @example DEC
Fokker_Planck = @decapode begin
  (ρ,Ψ)::Form0
  β⁻¹::Constant
  ∂ₜ(ρ) == ∘(⋆,d,⋆)(d(Ψ)∧ρ) + β⁻¹*Δ(ρ)
end
```

Specify the domain
``` @example DEC
spheremesh = loadmesh(Icosphere(6))
dualmesh = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(spheremesh);
subdivide_duals!(dualmesh, Barycenter())
```

Compile the simulation
``` @example DEC
simulation = eval(gensim(Fokker_Planck))
f = simulation(dualmesh, nothing)
```

Specify initial conditions. Ψ must be a smooth function. Choose an interesting eigenfunction. We require that ρ integrated over the surface is 1, since it is a PDF. On a sphere where ρ(x,y,z) is proportional to the x-coordinate, that means divide by 2π.
``` @example DEC
Δ0 = Δ(0, dualmesh)
Ψ = real.(eigs(Δ0, nev=32, which=:LR)[2][:,32])
ρ = map(point(dualmesh)) do (x,y,z)
  abs(x)
end / 2π
```

Let's define the structures which hold the constants and state variables for the
simulation, respectively.
``` @example DEC
constants_and_parameters = (β⁻¹ = 1e-2,)
u0 = ComponentArray(Ψ=Ψ, ρ=ρ)
```

Run the simulation.
``` @example DEC
tₑ= 20.0
problem = ODEProblem(f, u0, (0, tₑ), constants_and_parameters);
solution = solve(problem, Tsit5(), progress=true, progress_steps=1);
```

Verify that the probability distribbution function is still a probability distribution. We'll show that the sum of the values on the
dual 2-form integrate (sum to) unity,
``` @example DEC
s0 = dec_hodge_star(0, dualmesh)
@info sum(s0 * solution(tₑ).ρ)
@info any(solution(tₑ).ρ .≤ 0)
```

Now we will create a GIF.

``` @example DEC
function save_gif(file_name, soln)
  time = Observable(0.0)
  fig = Figure()
  Label(fig[1, 1, Top()], @lift("ρ at $($time)"), padding = (0, 0, 5, 0))
  ax = LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.mesh!(ax, spheremesh,
    color=@lift(soln($time).ρ),
    colorrange=(0,1),
    colormap=:jet)

  Colorbar(fig[1,2], msh)
  frames = range(0.0, tₑ; length=21)
  record(fig, file_name, frames; framerate = 10) do t
    time[] = t
  end
end
gif = save_gif("fokker_planck.gif", solution)
```

!["FokkerPlanck"](fokker_planck.gif)

