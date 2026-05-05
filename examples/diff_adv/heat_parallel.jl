using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using ComponentArrays
using OrdinaryDiffEq

Point3D = Point3{Float64}

rect = triangulated_grid(100, 100, 1, 1, Point3D)
d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(rect)
subdivide_duals!(d_rect, Circumcenter())

Heat = @decapode begin
  U::Form0
  ∂ₜ(U) == 100 * Δ(U)
end

simulate = evalsim(Heat)
fₘ = simulate(d_rect, nothing)

U = map(d_rect[:point]) do (x, _)
  x
end

u₀ = ComponentArray(U = U)

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4))
soln = solve(prob, KuttaPRK2p5(); dt = 1e-2)
soln.retcode != :Unstable || error("Solver was not stable")

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ))
soln = solve(prob, KuttaPRK2p5(); dt = 1e-2)
@info("Done")
