module Environment

using DiagrammaticEquations
using DiagrammaticEquations.Deca
using ..Canon
using Markdown

@docapode("Glens Law"
  ,"https://www.google.com"
  ,"Nye, J. F. (1957). The Distribution of Stress and Velocity in Glaciers and Ice-Sheets. Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 239(1216), 113–133. http://www.jstor.org/stable/100184"
  ,glen
  ,begin 
    Γ::Form1
    (A,ρ,g,n)::Constant

    Γ == (2/(n+2))*A*(ρ*g)^n
  end
)

@docapode("Halfar (Eq. 2)"
  ,"https://www.google.com"
  ,"Halfar, P. (1981), On the dynamics of the ice sheets, J. Geophys. Res., 86(C11), 11065–11072, doi:10.1029/JC086iC11p11065"
  ,halfar_eq2
  ,begin
    h::Form0
    Γ::Form1
    n::Constant

    ∂ₜ(h) == ∘(⋆, d, ⋆)(Γ  * d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end)


@docapode("Energy balance"
          ,"https://google.com"
          ,"energy balance equation from Budyko Sellers"
          ,energy_balance
          ,begin
  (Tₛ, ASR, OLR, HT)::Form0
  (C)::Constant

  Tₛ̇ == ∂ₜ(Tₛ) 

  Tₛ̇ == (ASR - OLR + HT) ./ C
end)

@docapode("Insolation"
          ,"https://google.com"
          ,""
          ,insolation
          , begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end)

@docapode("Warming"
          ,"https://google.com"
          ,""
          ,warming
          , begin
  (Tₛ)::Form0
  (A)::Form1

  A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.65e7)
end)

@docapode("Tracer"
          ,"https://google.com"
          ,""
          ,tracer
          , begin
  (c,C,F,c_up)::Form0
  (v,V,q)::Form1

  c_up == -1*⋆(L(v,⋆(c))) - ⋆(L(V,⋆(c))) - ⋆(L(v,⋆(C))) - ∘(⋆,d,⋆)(q) + F
end)

@docapode("Equation of State"
          ,"https://google.com"
          ,""
          ,equation_of_state
          , begin
  (b,T,S)::Form0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end)

@docapode("Boundary Conditions"
          ,"https://google.com"
          ,""
          ,boundary_conditions
          , begin
  (S,T)::Form0
  (Ṡ,T_up)::Form0
  v::Form1
  v_up::Form1
  Ṫ == ∂ₜ(T)
  Ṡ == ∂ₜ(S)
  v̇ == ∂ₜ(v)

  Ṫ == ∂_spatial(T_up)
  v̇ == ∂_noslip(v_up)
end)

end
