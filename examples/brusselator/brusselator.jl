using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
#using GLMakie
using Logging

using GeometryBasics: Point3
Point3D = Point3{Float64}

Brusselator = SummationDecapode(parse_decapode(
quote
  # Values living on vertices.
  (U, V, U2V, F, U̇, V̇)::Form0{X}
  # Scalars.
  (one, fourfour, threefour, α)::Parameter{X}
  # An intermediate variable.
  U2V == (U .* U) .* V
  # Specify how to compute the tangent variables.
  U̇ == one + U2V - (fourfour * U) + (α * Δ(U)) + F
  V̇ == (threefour * U) - U2V + (α * Δ(U))
  # Associate tangent variables with a state variable.
  U̇ == ∂ₜ(U)
  V̇ == ∂ₜ(V)
end))

# We resolve types of intermediate variables using sets of rules.
bespoke_op1_inf_rules = [
  (src_type = :Form0, tgt_type = :infer, replacement_type = :Form0, op = :Δ)]

bespoke_op2_inf_rules = [
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :.*),
  (proj1_type = :Form0, proj2_type = :Parameter, res_type = :infer, replacement_type = :Form0, op = :*),
  (proj1_type = :Parameter, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :*)]

infer_types!(Brusselator,
    vcat(bespoke_op1_inf_rules, op1_inf_rules_2D),
    vcat(bespoke_op2_inf_rules, op2_inf_rules_2D))

# Resolve overloads. i.e. map symbols to functions i.e. ~dispatch
resolve_overloads!(Brusselator)

# TODO: Create square domain of approximately 32x32 vertices.

# TODO: Create initial data.

# TODO: Create problem and run sim for t ∈ [0,1).

# TODO: Create problem and run sim for t ∈ [1.1,11.5].

