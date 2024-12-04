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

## Vorticity Formulation

```julia
eq11_inviscid_vorticity = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1

  𝐮 == d₁⁻¹(d𝐮)

  ∂ₜ(d𝐮) ==  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
```

## Initial Conditions

You can simulate this on either of the following initial conditions that are known to have periodic solutions for this inviscid formulation of NS.

### Taylor Vortices

![Plot of Taylor Vortex initial conditions]()

### Point Vortices

Our initial conditions here are Point vortices:

```julia
function point_vortex(pnt::Point3D, cntr::Point3D, p::PointVortexParams)
  gcd = great_circle_dist(pnt,cntr)
  p.τ / (cosh(3gcd/p.a)^2)
end
```

![Plot of Point Vortex initial conditions]()

## Expected Solutions

![Solution to Taylor Vortices GIF]()
![Solution to Point Vortices GIF]()

Based on the [configuration](config.toml), you can see different results that match the expected solutions from the literature.

Here is one set of results from using the inviscid Poisson formulation:

![Vorticity](vort.gif)

### Actual Solutions

These should be the unstable ones that come from the poor formulation.

![Solution to Taylor Vortices GIF]()
![Solution to Point Vortices GIF]()

## Laplacian Solver Fix

There are cohomological reasons why the above model formulation produces low-quality simualtions. The variable **X** is physically required to be in the kernel of $\Delta $, but that isn't guaranteed by the model formulation above. To fix this, you can use the **T** technique that requires a Laplacian solve as part of the update law.

This transformation can be implemented by editing the Decapode formulation and regenerating the simulator.

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


![Solution to Taylor Vortices GIF]()
![Solution to Point Vortices GIF]()


## Phenominological Assessment

These scenarios are used to test that a simulator achieves the correct phenomenology. In this case, we are looking for periodicity in the solution for vorticity. As the vortices advect around the sphere, there return to their original locations. This can be seen on the azimuth profile. The original formulation does not exhibit this phenomenon, but the corrected formulation does.

We can visualize the distribution of vorticity at the $\theta = 0.4$ latitude. The difference between the distributions at $t=0$ and $t=12$ is accumulated error.

![Azimuth Profile BAD](azimuth.png)
![Azimuth Profile GOOD](azimuth.png)


## Aritficial Viscosity

These simulations are intended to model an inviscid flow. To improve numerical stability, you can add a term for viscosity with a small, aphysical viscosity parameter just to allow error in the computed velocity profile to diffuse rather than accumulate.

```julia
eq11_inviscid_poisson = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1
  ψ::Form0

  ψ == Δ⁻¹(⋆(d𝐮))
  𝐮 == ⋆(d(ψ))

  ∂ₜ(d𝐮) ==  μ * ∘(⋆, d, ⋆, d)(d𝐮) +  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
```

```@example INFO
DocInfo.get_report(info) # hide
```
