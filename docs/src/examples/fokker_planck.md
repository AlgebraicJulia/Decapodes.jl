```@meta
EditURL = "../../literate/fokker_planck.jl"
```

We use here the formulation studied by Jordan, Kinderlehrer, and Otto in "The
Variational Formulation of the Fokker-Planck Equation" (1996).

The formulation they studied is that where the drift coefficient is the
gradient of (some potential) Ψ.

Load libraries.

````@example fokker_planck
using Catlab, CombinatorialSpaces, Decapodes, DiagrammaticEquations
using CairoMakie, ComponentArrays, LinearAlgebra, MLStyle, ComponentArrays, OrdinaryDiffEq
using GeometryBasics: Point3
Point3D = Point3{Float64}
using Arpack
````

Specify physics.

````@example fokker_planck
Fokker_Planck = @decapode begin
  (ρ,Ψ)::Form0
  β⁻¹::Constant
  ∂ₜ(ρ) == ∘(⋆,d,⋆)(d(Ψ)∧ρ) + β⁻¹*Δ(ρ)
end
````

Specify domain.

````@example fokker_planck
s = loadmesh(Icosphere(6));
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Barycenter());
nothing #hide
````

Compile.

````@example fokker_planck
sim = eval(gensim(Fokker_Planck))
fₘ = sim(sd, nothing)
````

Specify initial conditions.
Ψ must be a smooth function. Choose an interesting eigenfunction.

````@example fokker_planck
Δ0 = Δ(0,sd)
Ψ = real.(eigs(Δ0, nev=32, which=:LR)[2][:,32])
````

We require that ρ integrated over the surface is 1, since it is a PDF.
On a sphere where ρ(x,y,z) is proportional to the x-coordinate, that means divide by 2π.

````@example fokker_planck
ρ = map(point(sd)) do (x,y,z)
  abs(x)
end / 2π

constants_and_parameters = (β⁻¹ = 1e-2,)
u₀ = ComponentArray(Ψ=Ψ, ρ=ρ)
````

Run.

````@example fokker_planck
tₑ = 20.0
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters);
soln = solve(prob, Tsit5(), progress=true, progress_steps=1);
nothing #hide
````

Verify that the PDF is still a PDF.

````@example fokker_planck
s0 = dec_hodge_star(0,sd);
@info sum(s0 * soln(0).ρ)
@info sum(s0 * soln(tₑ).ρ) # ρ integrates to 1
@info any(soln(tₑ).ρ .≤ 0) # ρ is nonzero
````

Create GIF

````@example fokker_planck
function save_gif(file_name, soln)
  time = Observable(0.0)
  fig = Figure()
  Label(fig[1, 1, Top()], @lift("ρ at $($time)"), padding = (0, 0, 5, 0))
  ax = LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.mesh!(ax, s,
    color=@lift(soln($time).ρ),
    colorrange=(0,1),
    colormap=:jet)

  Colorbar(fig[1,2], msh)
  frames = range(0.0, tₑ; length=21)
  record(fig, file_name, frames; framerate = 10) do t
    time[] = t
  end
end
gif = save_gif("fokker_planck.gif", soln)

gif
````

using DisplayAs

DisplayAs.Text(DisplayAs.GIF(gif))

