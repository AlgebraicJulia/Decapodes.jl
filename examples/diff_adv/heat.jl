using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
Point2D = Point2{Float64}
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

U = map(d_rect[:point]) do (x,_)
        return x
    end

u₀ = ComponentArray(U=U)

constants_and_parameters = ()

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")
