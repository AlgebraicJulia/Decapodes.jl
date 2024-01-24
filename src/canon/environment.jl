using DiagrammaticEquations
using DiagrammaticEquations.Deca

@docapode(GlensLaw
  ,"link"
  ,"desc"
  ,glen
  ,begin 
    Γ::Form1
    (A,ρ,g,n)::Constant

    Γ == (2/(n+2))*A*(ρ*g)^n
  end
)

@docapode(halfar_eq2
  ,"link"
  ,"desc"
  ,halfar_eq2
  ,begin
    h::Form0
    Γ::Form1
    n::Constant

    ∂ₜ(h) == ∘(⋆, d, ⋆)(Γ * d(h) ∧ mag(♯(d(h)))^(n-1) ∧ (h^(n+2)))
  end
)

@docapode(KlausmeierEq2a
  ,"link"
  ,"desc"
  ,klausmeier_2a
  ,begin
    (n,w)::DualForm0
    dX::Form1
    (a,ν)::Constant
    ∂ₜ(w) == a - w - w * n^2 + ν * ℒ(dX, w)
  end
)

@docapode(KlausmeierEq2b
  ,"link"
  ,"desc"
  ,klausmeier_2b
  ,begin
    (n,w)::DualForm0
    m::Constant
    ∂ₜ(n) == w * n^2 - m*n + Δ(n)
  end)

@docapode(Lejeune
  ,"link"
  ,"desc"
  ,lejeune
  ,begin
    ρ::Form0
    (μ, Λ, L)::Constant

    ∂ₜ(ρ) == ρ * (1 - μ + (Λ - 1) * ρ - ρ*ρ) +
      0.5 * (L*L - ρ) * Δ(ρ) - 0.125 * ρ * Δ(ρ) * Δ(ρ)
  end)
