using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using LazyArrays: @~
using OrdinaryDiffEq

import Decapodes: default_dec_matrix_generate

# =============================================================================
# LazyArrays heat-equation proof of concept
#
# The gensim output for ∂ₜ(U) == κ * Δ(U) has two distinct phases:
#
#   SETUP (once, when simulate(mesh, operators) is called):
#     1. Fetch individual operator matrices:
#          ⋆₁, ⋆₀⁻¹, dual_d₁, d₀  via default_dec_matrix_generate
#     2. Contract them (sparse matrix-matrix products):
#          GenSim-M_GenSim-ConMat_0 = ⋆₀⁻¹ * dual_d₁ * ⋆₁ * d₀
#
#   PER-STEP (every ODE timestep):
#          mul!(•2, GenSim-M_GenSim-ConMat_0, U)   # one sparse mat-vec
#          U̇ .= κ .* •2
#
# The DEFAULT path materialises the Laplacian matrix once at setup time
# (sparse mat-mat products), then applies it cheaply at every timestep.
#
# The LAZYARRAYS path wraps each individual matrix in `@~` so the contraction
# step builds a lazy product instead of a dense/sparse materialisation.  At
# each timestep the lazy product applies the four operators in sequence, which
# avoids allocating the fully-contracted matrix but performs slightly more work
# per step.
# =============================================================================

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

# --- Default path: contracts operators into one matrix at setup time ----------
simulate_default = evalsim(Heat)
@info "=== DEFAULT: setup (contracts ⋆₀⁻¹ * dual_d₁ * ⋆₁ * d₀ once) ==="
@time fₘ_default = simulate_default(d_rect, nothing)

# --- LazyArrays path: keeps each operator matrix lazy; no contraction ---------
function lazy_generate(sd, my_symbol; hodge=GeometricHodge())
    M, op = default_dec_matrix_generate(sd, my_symbol, hodge)
    return (@~ M), op
end

simulate_lazy = evalsim(expand_operators(Heat))
@info "=== LAZYARRAYS: setup (no matrix-matrix products) ==="
@time fₘ_lazy = simulate_lazy(d_rect, lazy_generate)

# --- Initial conditions -------------------------------------------------------
U_initial = map(d_rect[:point]) do (x, _)
    x
end
u₀ = ComponentArray(U=U_initial)
constants_and_parameters = (κ=100.0,)

# Keep the same horizon as examples/diff_adv/heat.jl for direct comparison.
final_time = 11.5

# --- Solve: default -----------------------------------------------------------
prob_default = ODEProblem(fₘ_default, u₀, (0.0, final_time), constants_and_parameters)
@info "Precompiling default solver…"
solve(ODEProblem(fₘ_default, u₀, (0.0, 1e-4), constants_and_parameters), Tsit5())
@info "=== DEFAULT: solve ==="
@time soln_default = solve(prob_default, Tsit5())
soln_default.retcode != :Unstable || error("Default solver was not stable")

# --- Solve: LazyArrays --------------------------------------------------------
prob_lazy = ODEProblem(fₘ_lazy, u₀, (0.0, final_time), constants_and_parameters)
@info "Precompiling LazyArrays solver…"
solve(ODEProblem(fₘ_lazy, u₀, (0.0, 1e-4), constants_and_parameters), Tsit5())
@info "=== LAZYARRAYS: solve ==="
@time soln_lazy = solve(prob_lazy, Tsit5())
soln_lazy.retcode != :Unstable || error("LazyArrays solver was not stable")
