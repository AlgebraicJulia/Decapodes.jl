using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using LazyArrays: @~
using OrdinaryDiffEq

import Decapodes: default_dec_matrix_generate

Point3D = Point3{Float64}

rect = triangulated_grid(100, 100, 1, 1, Point3D)
d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(rect)
subdivide_duals!(d_rect, Circumcenter())

Heat = @decapode begin
    U::Form0
    κ::Constant
    # Δ(U) expands to ⋆₀⁻¹ ∘ dual_d₁ ∘ ⋆₁ ∘ d₀.
    ∂ₜ(U) == κ * Δ(U)
end

function lazy_generate(sd, my_symbol; hodge=GeometricHodge())
    M, op = default_dec_matrix_generate(sd, my_symbol, hodge)
    return (@~ M), op
end

simulate = evalsim(expand_operators(Heat))
fₘ = simulate(d_rect, lazy_generate)

U_initial = map(d_rect[:point]) do (x_coord, _y_coord)
    x_coord
end
u₀ = ComponentArray(U=U_initial)
constants_and_parameters = (κ=100.0,)

final_time = 11.5
prob = ODEProblem(fₘ, u₀, (0.0, final_time), constants_and_parameters)
soln = solve(prob, Tsit5())

soln.retcode != :Unstable || error("Solver was not stable")
