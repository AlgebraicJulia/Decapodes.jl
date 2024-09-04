# Harmonics of the Sphere

This page shows how to use Decapodes tooling to explore the harmonics of a discrete manifold. This isn't using any decapode specific code, but it 
is emblematic of a more advanced analysis you might want to do on your decapode.

In this case we are trying to visualize the roots of the Laplacian on a discrete manifold.

Load the dependencies

```@example Harmonics
# Meshing:
using CombinatorialSpaces
using CoordRefSystems
using GeometryBasics: Point3
const Point3D = Point3{Float64};

# Visualization:
using CairoMakie

# Simulation:
using LinearAlgebra
```

Load the mesh

```@example Harmonics
const RADIUS = 1.0
s = loadmesh(Icosphere(3, RADIUS));
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Barycenter());
```

Compute the laplacian eigenvectors using [LinearAlgebra.eigen](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.eigen). This requires making the sparse Laplacian matrix dense with `collect`. Alternatively, use [Arpack.jl](https://arpack.julialinearalgebra.org/stable/).

```@example Harmonics
Δ0 = -Δ(0,sd)
λ = eigen(collect(Δ0))
```

Let's check that our eigenvalues satisfy the right equation.
The first eigenvector should be the kernel of the laplacian.
So the following norm should be close to 0.

```@example Harmonics
q1 = λ.vectors[:,1]
norm(Δ0 *q1)
```

The first eigenvector is boring to visualize, because it is constant.
So we will make some plots of the second eigenvector.
If you run this on the desktop, you can use GLMakie and get an interactive plot to explore. We will just draw two angles.

```@example Harmonics
q = λ.vectors[:,2]
fig = Figure()
Label(fig[1, 1, Top()], "Default Angle", padding = (0, 0, 5, 0))
ax = LScene(fig[1,1], scenekw=(lights=[],))
msh = CairoMakie.mesh!(ax, s, color=q)
Colorbar(fig[1,2], msh, size=32)
# Second Angle
Label(fig[2, 1, Top()], "Bottom Angle", padding = (0, 0, 5, 0))
ax = LScene(fig[2,1], scenekw=(lights=[],))
update_cam!(ax.scene, Vec3f(-1/2,-1/2,1.0/2), Vec3f(1,1,1), Vec3f(0, 0, 1))
msh = CairoMakie.mesh!(ax, s, color=q)
Colorbar(fig[2,2], msh, size=32)
fig
```

```@example Harmonics
q = λ.vectors[:,12]
fig = Figure()
Label(fig[1, 1, Top()], "Default Angle", padding = (0, 0, 5, 0))
ax = LScene(fig[1,1], scenekw=(lights=[],))
msh = CairoMakie.mesh!(ax, s, color=q)
Colorbar(fig[1,2], msh, size=32)
# Second Angle
Label(fig[2, 1, Top()], "Bottom Angle", padding = (0, 0, 5, 0))
ax = LScene(fig[2,1], scenekw=(lights=[],))
update_cam!(ax.scene, Vec3f(-1/2,-1/2,1.0/2), Vec3f(1,1,1), Vec3f(0, 0, 1))
msh = CairoMakie.mesh!(ax, s, color=q)
Colorbar(fig[2,2], msh, size=32)
fig
```

# Exploring solutions with Krylov methods

We can also use the information about the eigenvectors for spectral techniques in solving the equations. Krylov methods are a bridge between
linear solvers and spectral information.

```julia
using Krylov
b = zeros(nv(sd))
b[1] = 1
b[end] = -1
x, stats = Krylov.gmres(Δ0, b, randn(nv(sd)), restart=true, memory=20, atol = 1e-10, rtol=1e-8, history=true, itmax=10000)
x̂ = x .- sum(x)./length(x)
norm(x̂)
stats
norm(Δ0*(x) - b)
```
