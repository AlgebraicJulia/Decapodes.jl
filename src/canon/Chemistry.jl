module Chemistry

using DiagrammaticEquations
using DiagrammaticEquations.Deca
using ..Canon
using Markdown

@docapode("Brusselator"
  ,"https://en.wikipedia.org/wiki/Brusselator"
  ,"A model of reaction-diffusion for an oscillatory chemical system."
  ,brusselator
  ,begin
    # Values living on vertices.
    (U, V)::Form0 # State variables.
    U2V::Form0 # Named intermediate variables.
    (U̇, V̇)::Form0 # Tangent variables.
    # Scalars.
    (α)::Constant
    F::Parameter
    # A named intermediate variable.
    U2V == (U .* U) .* V
    # Specify how to compute the tangent variables.
    U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
    V̇ == (3.4 * U) - U2V + (α * Δ(V))
    # Associate tangent variables with a state variable.
    ∂ₜ(U) == U̇
    ∂ₜ(V) == V̇
end)

@docapode("Gray-Scott"
  ,""
  ,"A model of reaction-diffusion"
  ,GrayScott
  ,begin
    (U, V)::Form0
    (UV2)::Form0
    (U̇, V̇)::Form0
    (f, k, rᵤ, rᵥ)::Constant

    UV2 == (U .* (V .* V))
    U̇ == rᵤ * Δ(U) - UV2 + f * (1 .- U)
    V̇ == rᵥ * Δ(V) + UV2 - (f + k) .* V
    ∂ₜ(U) == U̇
    ∂ₜ(V) == V̇
  end
)

end
