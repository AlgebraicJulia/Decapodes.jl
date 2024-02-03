using DiagrammaticEquations
using DiagrammaticEquations.Deca

@docapode(GlensLaw
  ,"https://www.google.com"
  ,"Nye, J. F. (1957). The Distribution of Stress and Velocity in Glaciers and Ice-Sheets. Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 239(1216), 113–133. http://www.jstor.org/stable/100184"
  ,glen
  ,begin 
    Γ::Form1
    (A,ρ,g,n)::Constant

    Γ == (2/(n+2))*A*(ρ*g)^n
  end
)

@docapode(halfar_eq2
  ,"https://www.google.com"
  ,"Halfar, P. (1981), On the dynamics of the ice sheets, J. Geophys. Res., 86(C11), 11065–11072, doi:10.1029/JC086iC11p11065"
  ,halfar_eq2
  ,beginhalfar_eq2 = @decapode begin
    h::Form0
    Γ::Form1
    n::Constant

    ∂ₜ(h) == ∘(⋆, d, ⋆)(Γ  * d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end

  end
)

@docapode(KlausmeierEq2a
  ,"https://www.google.com"
  ,"Klausmeier, CA. “Regular and irregular patterns in semiarid vegetation.” Science (New York, N.Y.) vol. 284,5421 (1999): 1826-8. doi:10.1126/science.284.5421.1826"
  ,klausmeier_2a
  ,begin
    (n,w)::DualForm0
    dX::Form1
    (a,ν)::Constant
    ∂ₜ(w) == a - w - w * n^2 + ν * ℒ(dX, w)
  end
)

@docapode(KlausmeierEq2b
  ,"https://www.google.com"
  ,"...

  Klausmeier, CA. “Regular and irregular patterns in semiarid vegetation.” Science (New York, N.Y.) vol. 284,5421 (1999): 1826-8. doi:10.1126/science.284.5421.1826"
  ,klausmeier_2b
  ,begin
    (n,w)::DualForm0
    m::Constant
    ∂ₜ(n) == w * n^2 - m*n + Δ(n)
  end)

@docapode(Lejeune
  ,"https://www.google.com"
  ,"Lejeune, O., & Tlidi, M. (1999). A Model for the Explanation of Vegetation Stripes (Tiger Bush). Journal of Vegetation Science, 10(2), 201–208. https://doi.org/10.2307/3237141"
  ,lejeune
  ,begin
    ρ::Form0
    (μ, Λ, L)::Constant

    ∂ₜ(ρ) == ρ * (1 - μ + (Λ - 1) * ρ - ρ*ρ) +
      0.5 * (L*L - ρ) * Δ(ρ) - 0.125 * ρ * Δ(ρ) * Δ(ρ)
  end)
