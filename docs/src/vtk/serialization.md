# Simulation Data Serialization

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..", "docinfo.jl"))
info = DocInfo.Info()
```

Decapodes simulation results can be serialized to the [VTK](https://vtk.org/) format for visualization in [ParaView](https://www.paraview.org/), the industry-standard tool for scientific visualization. This page demonstrates how to export simulation data using [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).

## Dependencies

```@example DEC
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis

using Catlab
using CombinatorialSpaces
using ComponentArrays
using DiagrammaticEquations
using Decapodes
using Distributions
using LinearAlgebra
using OrdinaryDiffEq
using WriteVTK

using GeometryBasics: Point3
Point3D = Point3{Float64}
nothing # hide
```

## Define the Model

We define a simple diffusion model to demonstrate serialization. This models the spread of a concentration `C` over time, governed by Fick's law.

```@example DEC
Diffusion = @decapode begin
  (C, Ċ)::Form0
  ϕ::Form1
  k::Constant

  # Fick's first law
  ϕ == k * (d₀(C))

  # Diffusion equation
  Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
  ∂ₜ(C) == Ċ
end

to_graphviz(Diffusion)
```

## The Mesh

We create a `triangulated_grid` mesh and its dual complex. The primal mesh `s` stores the geometry (vertices, edges, triangles), while the dual complex `sd` is used for DEC operators.

```@example DEC
s = triangulated_grid(30, 10, 1, 1, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
wireframe!(ax, s)
save("serialization_mesh.png", fig)
nothing # hide
```

!["Triangulated Grid Mesh"](serialization_mesh.png)

## Initial Conditions and Simulation

We set up a Gaussian initial condition and solve the diffusion equation.

```@example DEC
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in s[:point]]

fig = Figure()
ax = CairoMakie.Axis(fig[1,1], aspect = AxisAspect(3.0))
msh = mesh!(ax, s; color=c, colorrange=extrema(c))
Colorbar(fig[1,2], msh)
save("serialization_initial.png", fig)
nothing # hide
```

!["Initial Conditions"](serialization_initial.png)

```@example DEC
sim = eval(gensim(Diffusion))
fₘ = sim(sd, nothing, DiagonalHodge())

u₀ = ComponentArray(C=c)
constants = ComponentArray(k=0.05)

tₑ = 100
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), constants)
soln = solve(prob, Tsit5())
soln.retcode
```

## Exporting to VTK

### Mesh Preparation

To write VTK files, we first extract the mesh geometry into the format expected by WriteVTK.jl. The points are stored as a `3 × N` matrix, and each triangle is represented as a `MeshCell`.

```@example DEC
points = hcat([collect(p) for p in s[:point]]...)
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, triangle_vertices(s, t))
         for t in 1:ntriangles(s)]
nothing # hide
```

### Saving a 0-form (Scalar Point Data)

A single snapshot of the simulation can be saved as a `.vtu` (VTK Unstructured Grid) file. A 0-form like `C` has one value per vertex, so it is attached as point data with `VTKPointData`.

```@example DEC
vtk = vtk_grid("diffusion_final", points, cells)
vtk["C", VTKPointData()] = soln(tₑ).C
vtk_save(vtk)
```

### Saving a 1-form (Vector Point Data)

A 1-form has one scalar value per edge. To visualize it in ParaView, we use the **sharp** operator (`♯`) to convert it into a vector field at vertices. The `AltPPSharp` interpolation method produces a primal-primal sharp suitable for visualization on 2D simplicial complexes.

We first compute the flux 1-form `ϕ = k * d₀(C)` from the solution, then sharpen it. The resulting vector field has one vector per vertex and is exported as a `3 × N` matrix with `VTKPointData`. We define a small helper `vec_to_mat` to convert the vector of points into the matrix layout that WriteVTK expects.

```@example DEC
vec_to_mat(vecs) = hcat(collect.(vecs)...)

d₀_mat = dec_differential(0, sd)
ϕ = 0.05 * d₀_mat * soln(tₑ).C

sharp = ♯_mat(sd, AltPPSharp())
ϕ_vec = sharp * ϕ

vtk = vtk_grid("diffusion_flux", points, cells)
vtk["C", VTKPointData()] = soln(tₑ).C
vtk["flux", VTKPointData()] = vec_to_mat(ϕ_vec)
vtk_save(vtk)
```

### Saving a 2-form (Cell Data)

A 2-form has one scalar value per triangle (2-simplex). In VTK terminology, triangles are **cells**, so 2-form data is attached with `VTKCellData`.

Here we compute a 2-form by taking the exterior derivative `d₁` of the flux 1-form. On a flat mesh, `d₁ ∘ d₀ = 0` exactly, so these values will be near-zero — this serves as a numerical consistency check. We also compute a more visual example: the average vertex concentration per triangle, which illustrates how any per-triangle quantity is exported as cell data.

```@example DEC
d₁_mat = dec_differential(1, sd)
curl_ϕ = d₁_mat * ϕ

C_final = soln(tₑ).C
avg_C = [sum(C_final[triangle_vertices(s, t)]) / 3
         for t in 1:ntriangles(s)]

vtk = vtk_grid("diffusion_celldata", points, cells)
vtk["C", VTKPointData()] = C_final
vtk["curl_flux", VTKCellData()] = curl_ϕ
vtk["avg_C", VTKCellData()] = avg_C
vtk_save(vtk)
```

### Time Series with ParaView Collection

For visualizing the full time evolution in ParaView, we export each timestep as a separate `.vtu` file and collect them into a `.pvd` (ParaView Data) file. ParaView can then load the `.pvd` file and animate through the timesteps. This example exports the 0-form, the sharpened 1-form, and the 2-form together.

```@example DEC
pvd = paraview_collection("diffusion_series")
timestamps = range(0.0, tₑ, length=50)
for (i, t) in enumerate(timestamps)
  local C_t = soln(t).C
  local ϕ_t = 0.05 * d₀_mat * C_t
  local avg_C_t = [sum(C_t[triangle_vertices(s, t_idx)]) / 3
                   for t_idx in 1:ntriangles(s)]

  local vtk_step = vtk_grid("diffusion_series_$i", points, cells)
  vtk_step["C", VTKPointData()] = C_t
  vtk_step["flux", VTKPointData()] = vec_to_mat(sharp * ϕ_t)
  vtk_step["avg_C", VTKCellData()] = avg_C_t
  collection_add_timestep(pvd, vtk_step, t)
end
vtk_save(pvd)
```

The resulting `diffusion_series.pvd` file can be opened in ParaView to visualize the diffusion process over time with full playback controls.

## Cleanup

```@example DEC
foreach(filter(endswith(".vtu"), readdir())) do f # hide
  rm(f) # hide
end # hide
foreach(filter(endswith(".pvd"), readdir())) do f # hide
  rm(f) # hide
end # hide
nothing # hide
```

```@example INFO
DocInfo.get_report(info) # hide
```
