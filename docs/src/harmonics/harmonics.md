# Harmonics of the Sphere

This page shows how to use Decapodes tooling to explore the harmonics of a discrete manifold. This isn't using any decapode specific code, but it 
is emblematic of a more advanced analysis you might want to do on your decapode.

In this case we are trying to visualize the roots of the Laplacian on a discrete manifold.

Load the dependencies

```@example Harmonics
using CombinatorialSpaces
# using Decapodes
# using DiagrammaticEquations
using CoordRefSystems
using GeometryBasics: Point3
Point3D = Point3{Float64};
# Visualization:
using CairoMakie
# Simulation:
# using ComponentArrays
using LinearAlgebra
# using LinearAlgebra: factorize
Point3D = Point3{Float64};
```

Load the mesh

```@example Harmonics
BARYCENTRIC = true
const RADIUS = 1.0
s = loadmesh(Icosphere(3, RADIUS));
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Barycenter());
```

Compute the laplacian eigenvectors

```@example Harmonics
Δ0 = -Δ(0,sd)
λ = eigen(collect(Δ0))
```

```@example Harmonics
λ.values
q1 = λ.vectors[:,1]
norm(Δ0 *q1)
```

Now make the plot

```@example Harmonics
q = λ.vectors[:,15]
fig = Figure()
Label(fig[1, 1, Top()], "eigenvector", padding = (0, 0, 5, 0))
ax = LScene(fig[1,1], scenekw=(lights=[],))
update_cam!(ax.scene, Vec3f(0,0,1.0), Vec3f(0,0,0), Vec3f(0, 1, 1))
msh = GLMakie.mesh!(ax, s, color=q)
Colorbar(fig[1,2], msh, size=32)
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
# x, stats = Krylov.gmres(Δ0, b, randn(nv(sd)), restart=true, memory=100, atol = 1e-10, rtol=1e-8)
# x, stats = Krylov.cg(Δ0, b, ones(nv(sd)), atol = 1e-10, rtol=1e-8, history=true)
# x, stats = Krylov.bicgstab(Δ0, b, randn(length(b)), restart=true, memory=45, atol = 1e-14, rtol=1e-8)
x̂ = x .- sum(x)./length(x)
norm(x̂)
stats
norm(Δ0*(x) - b)
```
