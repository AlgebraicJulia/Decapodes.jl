# FAQ

## 1. How do I incorporate scalar or vector field input data where you have a function of the embedded coordinates?

We can take a look at the [Brusselator page](../brussel/brussel.md#initial-data) which sets the values of each point on its mesh to a value as determined by some function. This can be done in a similar manner in 1D and 2D. 

The Brusselator also demonstrates, with the variable `F`, how one can go about changing the function by which these values are set over time.

## 2. How do I incorporate input data from a file with linear interpolation?

The Grigoriev Ice Cap model has a section where after the initial data is loaded from a TIF, the data is interpolated so that it may fit over a discrete mesh of our choosing. The link for that is [here](../grigoriev/grigoriev.md#loading-a-scientific-dataset).

## 3. How do I set boundary conditions like fixed value, no-flux, and no-slip?

Boundary conditions can be set by using "collages", which can take two variables among two different Decapodes and apply a function on the first. A general workflow would be to have the first Decapode encode the physics and have the second one encode values for boundary conditions. They can be related by a function that will mask the first variable and replace the desired values with second. An example of applying fixed boundary conditions would be in the [Brusselator page](../brussel/brussel.md#boundary-conditions).

A similar workflow can be used for "no-flux" and "no-slip" conditions by fixing the value of the appropriate variable to be 0. 

## 4. How do I plot derived quantities?

Plotting in DECAPODES is commonly done with the [Makie](https://github.com/MakieOrg/Makie.jl) package in Julia. Makie allows for creating both still images, which are useful for visualizing the mesh itself and initial/final conditions, and videos, which can capture the full simulation from start to end.

- For [1D visualization](../ice_dynamics/ice_dynamics.md#visualize)

- For [2D visualization](../ice_dynamics/ice_dynamics.md#visualize-2d)

- For [3D visualization](../ice_dynamics/ice_dynamics.md#2-manifold-in-3d)


## 5. How do I add artificial diffusion for 0- or 1-forms?

Without viscosity - i.e. when ``μ = 0`` - the incompressible (inviscid) Navier-Stokes equations can be formulated like so:

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

Adding a viscosity term can be accomplished by simply added the appropriate term, and declaring the ``μ`` constant:

```julia
eq11_viscid_poisson = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1
  ψ::Form0
  μ::Constant

  ψ == Δ⁻¹(⋆(d𝐮))
  𝐮 == ⋆(d(ψ))

  ∂ₜ(d𝐮) ==  μ * ∘(⋆, d, ⋆, d)(d𝐮) + (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
```

More demonstrations on how to iterate between formulations of the same physics (the incompressible Navier-Stokes equations) is available in further detail on the [Vortices](../navier_stokes/ns.md) docs page and in the script available there.

## 6. How do I use multigrid methods?

To use multigrid methods in the Laplacian solver, you need to create a `PrimalGeometricMapSeries` that will take a coarse mesh and apply a subdivision method to it some number of times. After that, just use this result as you would a regular mesh for simulation.

```julia
s = triangulated_grid(100,100,10,10,Point3D)

# Binary subdivide 4 times
series = PrimalGeometricMapSeries(s, binary_subdivision_map, 4);

# Retrieve highest resolution mesh
sd = finest_mesh(series)

  ...

f_mg = sim_mg(series, generate);
```

## 7. What are general workflows for DECAPODES?

A common workflow is to iterate through multiple different models as is done in the [Vorticity Model page](../navier_stokes/ns.md). A formulation is first done with a direct vorticity formulation but a quick run finds that this setup is unstable. A second formulation introduces a Laplacian solve which produces nice results.

Similar workflows may retain the same model but may iterate on the types of meshes/initial conditions used. An excellent example of this is found in the [Glacial Flow page](../ice_dynamics/ice_dynamics.md) where the model is first run in a [1D](../ice_dynamics/ice_dynamics.md#Define-a-mesh) setting and then quickly promoted to both [2D](../ice_dynamics/ice_dynamics.md#Define-our-mesh) and [3D](../ice_dynamics/ice_dynamics.md#2-Manifold-in-3D). This allows either running some dynamics in a more complicated setting, as just discussed, or allows for simplifying higher dimensional models by some sort of symmetry.
