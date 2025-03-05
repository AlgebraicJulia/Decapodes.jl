module Biology

using DiagrammaticEquations
using DiagrammaticEquations.Deca
using ..Canon
using Markdown

@docapode("Klausmeier (Eq. 2a)"
  ,""
  ,"Klausmeier, CA. “Regular and irregular patterns in semiarid vegetation.” Science (New York, N.Y.) vol. 284,5421 (1999): 1826-8. doi:10.1126/science.284.5421.1826"
  ,klausmeier_2a
  ,begin
    (n,w)::DualForm0
    dX::Form1
    (a,ν)::Constant
    ∂ₜ(w) == a - w - w * n^2 + ν * ℒ(dX, w)
  end
)

@docapode("Klausmeier (Eq. 2b)"
  ,""
  ,"ibid."
  ,klausmeier_2b
  ,begin
    (n,w)::DualForm0
    m::Constant
    ∂ₜ(n) == w * n^2 - m*n + Δ(n)
  end)

@docapode("Lejeune"
  ,""
  ,"Lejeune, O., & Tlidi, M. (1999). A Model for the Explanation of Vegetation Stripes (Tiger Bush). Journal of Vegetation Science, 10(2), 201–208. https://doi.org/10.2307/3237141"
  ,lejeune
  ,begin
    ρ::Form0
    (μ, Λ, L)::Constant

    ∂ₜ(ρ) == ρ * (1 - μ + (Λ - 1) * ρ - ρ*ρ) +
      0.5 * (L*L - ρ) * Δ(ρ) - 0.125 * ρ * Δ(ρ) * Δ(ρ)
  end)


@docapode("Kealy"
          ,""
          ,""
          ,kealy
          ,begin
            (n,w)::DualForm0
            dX::Form1
            (a,ν)::Constant
            ∂ₜ(w) == a - w - w * n^2 + ν * Δ(w)
  end
)

@docapode("Turing Continuous Ring"
  ,""
  ,""
  ,turing_continuous_ring
  ,begin
    (X, Y)::Form0
    (μ, ν, a, b, c, d)::Constant

    ∂ₜ(X) == a*X + b*Y + μ*Δ(X)
    ∂ₜ(Y) == c*X + d*Y + ν*Δ(X)
end)

end
