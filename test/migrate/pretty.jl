using Test
using MLStyle
using Decapodes
using Decapodes.decapodes

@testset "Pretty Printing" begin
mdl = parse_decapode(quote
  (Q, Tₛ, ASR, OLR, HT)::Form0
  (α, A,B,C)::Constant
  (D,cosϕᵖ,cosϕᵈ)::Constant

  Tₛ̇ == ∂ₜ(Tₛ) 
  Tₛ̇ == (ASR - OLR + HT) / C
  ASR == (1 - α) * Q
  OLR == A + (B * Tₛ)

  HT == (D ./ cosϕᵖ) * ⋆(d(cosϕᵈ * ⋆(d(Tₛ))))
end)
output = "Context:\n  Q::Form0 over I\n  Tₛ::Form0 over I\n  ASR::Form0 over I\n  OLR::Form0 over I\n  HT::Form0 over I\n  α::Constant over I\n  A::Constant over I\n  B::Constant over I\n  C::Constant over I\n  D::Constant over I\n  cosϕᵖ::Constant over I\n  cosϕᵈ::Constant over I\nEquations:\nTₛ̇   = ∂ₜ(Tₛ)\nTₛ̇   = ASR - OLR + HT / C\nASR   = 1 - α * Q\nOLR   = A + B * Tₛ\nHT   = D ./ cosϕᵖ * ⋆(d(cosϕᵈ * ⋆(d(Tₛ))))\n" 

# output = """Context:
#   Q::Form0 over I
#   Tₛ::Form0 over I
#   ASR::Form0 over I
#   OLR::Form0 over I
#   HT::Form0 over I
#   α::Constant over I
#   A::Constant over I
#   B::Constant over I
#   C::Constant over I
#   D::Constant over I
#   cosϕᵖ::Constant over I
#   cosϕᵈ::Constant over I
# Equations:
# Tₛ̇   = ∂ₜ(Tₛ)
# Tₛ̇   = ASR - OLR + HT / C
# ASR   = 1 - α * Q
# OLR   = A + B * Tₛ
# HT   = D ./ cosϕᵖ * ⋆(d(cosϕᵈ * ⋆(d(Tₛ))))"""

@test sprint((io, x)->Decapodes.pprint(io, x), mdl) == output
end