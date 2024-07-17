# Navier Stokes Vorticity Model

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..", "docinfo.jl"))
info = DocInfo.Info()
```

This is a discretization of the incompressible Navier Stokes equations using the Discrete Exterior Calculus.

The formulations are based on those given by [Mohamed, Hirani, Samtaney](https://arxiv.org/abs/1508.01166) (in turn from [Marsden, Ratiu, Abraham](https://link.springer.com/book/10.1007/978-1-4612-1029-0)).

However, different choices in discretization are chosen for purposes of brevity, to demonstrate novel discretizations of certain operators, and to demonstrate the automated Decapodes workflow.

The full code that generated these results is available in [a julia script](ns.jl).

We give the vorticity formulation of the inviscid incompressible Navier-Stokes momentum equation as follows:

```julia
eq11_inviscid_poisson = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1
  ψ::Form0

  ψ == Δ⁻¹(⋆(d𝐮))
  𝐮 == ⋆(d(ψ))

  ∂ₜ(d𝐮) ==  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
```

Our initial conditions here are Point vortices:

```julia
function point_vortex(pnt::Point3D, cntr::Point3D, p::PointVortexParams)
  gcd = great_circle_dist(pnt,cntr)
  p.τ / (cosh(3gcd/p.a)^2)
end
```

Based on the [configuration](config.toml), you can see different results that match the expected solutions from the literature.

Here is one set of results from using the inviscid Poisson formulation:

![Vorticity](vort.gif)

We can visualize the distribution of vorticity at the $\theta = 0.4$ latitude. The difference between the distributions at $t=0$ and $t=12$ is accumulated error.

![Azimuth Profile](azimuth.png)

```@example INFO
DocInfo.get_report(info) # hide
```
