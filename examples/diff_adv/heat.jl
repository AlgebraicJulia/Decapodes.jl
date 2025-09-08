using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using Distributions
using CairoMakie

lx = ly = 10
s = triangulated_grid(lx, ly, 0.1, 0.1, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

Heat = @decapode begin
    T::Form0
    D::Constant
    ∂ₜ(T) == D * Δ(T)
end

simulate = evalsim(Heat)
  
fₘ = simulate(sd, nothing)

T_dist = MvNormal([lx/2, ly/2], [1, 1])
T = [pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]

D = 0.1

u₀ = ComponentArray(T=T)

constants_and_parameters = (D = D,)

tₑ = 100

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=soln.u[end].T, colormap=:jet)
Colorbar(fig[1,2], msh)
save("2d_heat.png", fig)