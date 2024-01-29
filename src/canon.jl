module Canon

Brusselator_expr = quote
  # Values living on vertices.
  (U, V)::Form0{X} # State variables.
  (U2V, One)::Form0{X} # Named intermediate variables.
  (U̇, V̇)::Form0{X} # Tangent variables.
  # Scalars.
  (α)::Constant{X}
  F::Parameter{X}
  # A named intermediate variable.
  U2V == (U .* U) .* V
  # Specify how to compute the tangent variables.
  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(U))
  # Associate tangent variables with a state variable.
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

"""    Brusselator

The model of reaction diffusion for an oscillatory chemical system.

[Source](https://en.wikipedia.org/wiki/Brusselator)

Model:

$(Brusselator_expr)
"""
Brusselator = parse_decapode(Brusselator_expr)

GrayScott = @decapode begin
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

oscillator = @decapode begin
  X::Form0
  V::Form0
  k::Constant

  ∂ₜ(X) == V
  ∂ₜ(V) == -k*(X)
end

energy_balance = @decapode begin
  (Tₛ, ASR, OLR, HT)::Form0
  (C)::Constant

  Tₛ̇ == ∂ₜ(Tₛ) 

  Tₛ̇ == (ASR - OLR + HT) ./ C
end

absorbed_shortwave_radiation = @decapode begin
  (Q, ASR)::Form0
  α::Constant

  ASR == (1 .- α) .* Q
end

outgoing_longwave_radiation = @decapode begin
  (Tₛ, OLR)::Form0
  (A,B)::Constant

  OLR == A .+ (B .* Tₛ)
end

heat_transfer = @decapode begin
  (HT, Tₛ)::Form0
  (D,cosϕᵖ,cosϕᵈ)::Constant

  HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
end

insolation = @decapode begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end

warming = @decapode begin
  (Tₛ)::Form0
  (A)::Form1

  A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.65e7)
end

momentum = @decapode begin
  (f,b)::Form0
  (v,V,g,Fᵥ,uˢ,v_up)::Form1
  τ::Form2
  U::Parameter

  uˢ̇ == ∂ₜ(uˢ)

  v_up == -1 * L(v,v) - L(V,v) - L(v,V) -
       f∧v - ∘(⋆,d,⋆)(uˢ)∧v - d(p) + b∧g - ∘(⋆,d,⋆)(τ) + uˢ̇ + Fᵥ

  uˢ̇ == force(U)
end

tracer = @decapode begin
  (c,C,F,c_up)::Form0
  (v,V,q)::Form1

  c_up == -1*⋆(L(v,⋆(c))) - ⋆(L(V,⋆(c))) - ⋆(L(v,⋆(C))) - ∘(⋆,d,⋆)(q) + F
end

equation_of_state = @decapode begin
  (b,T,S)::Form0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end

boundary_conditions = @decapode begin
  (S,T)::Form0
  (Ṡ,T_up)::Form0
  v::Form1
  v_up::Form1
  Ṫ == ∂ₜ(T)
  Ṡ == ∂ₜ(S)
  v̇ == ∂ₜ(v)

  Ṫ == ∂_spatial(T_up)
  v̇ == ∂_noslip(v_up)
end

Ficks_Law = @decapode begin
  C::Form0
  ϕ::Form1

  # Fick's first law
  ϕ ==  k(d₀(C))
end

Advection = @decapode begin
  C::Form0
  (ϕ,V)::Form1

  ϕ == ∧₀₁(C,V)
end

Superposition = @decapode begin
  (C, Ċ)::Form0
  (ϕ, ϕ₁, ϕ₂)::Form1

  ϕ == ϕ₁ + ϕ₂
  Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
  ∂ₜ(C) == Ċ
end

Lie = @decapode begin
  C::Form0
  V::Form1
  dX::Form1

  V == ∘(⋆,⋆)(C ∧ dX)
end

Diffusion = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
  # Diffusion equation
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end

Diffusion = @decapode DiffusionQuantities begin
  C::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
end

Advection = @decapode DiffusionQuantities begin
  C::Form0{X}
  (V, ϕ)::Form1{X}
  ϕ == ∧₀₁{X}(C,V)
end

Superposition = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  (ϕ, ϕ₁, ϕ₂)::Form1{X}

  ϕ == ϕ₁ + ϕ₂
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end

Poiseuille = @decapode begin
  P::Form0
  q::Form1
  (R, μ̃ )::Constant

  # Laplacian of q for the viscous effect
  Δq == Δ(q)
  # Gradient of P for the pressure driving force
  ∇P == d(P)

  # Definition of the time derivative of q
  ∂ₜ(q) == q̇

  # The core equation
  q̇ == μ̃  * ∂q(Δq) + ∇P + R * q
end

Poiseuille_Density = @decapode begin
  q::Form1
  (P, ρ)::Form0
  (k, R, μ̃ )::Constant

  # Poiseuille Flow
  ∂ₜ(q) == q̇
  ∇P == d(P)
  q̇ == μ̃ * ∂q(Δ(q)) - ∇P + R * q
  
  # Pressure/Density Coupling
  P == k * ρ
  ∂ₜ(ρ) == ρ̇
  ρ_up == ∘(⋆, d, ⋆)(-1 * ∧₀₁(ρ,q)) # advection
  
  # Boundary conditions
  ρ̇ == ∂ρ(ρ_up)
end

halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ∂ₜ(h) == ∘(⋆, d, ⋆)(Γ * d(h) ∧ mag(♯(d(h)))^(n-1) ∧ (h^(n+2)))
end

glens_law = @decapode begin
  Γ::Form1
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

NavierStokes = @decapode begin
  (V, V̇, G)::Form1{X}
  (ρ, ṗ, p)::Form0{X}
  
  Δ₁(V) + third(d₀(δ₁(V)))
  
  ∂ₜ(V) == neg₁(L₁′(V, V)) + 
    kᵥ(Δ(V) + (1/3)*d(δ(V))) / avg₀₁(ρ) +
    d(0.5 * (i₁′(V, V))) +
    -1 * d(p) / avg₀₁(ρ) +
    G

  ∂ₜ(p) == -1 * ⋆(L₀(V, ⋆(p)))
end

Schroedinger = @decapode begin
  (i,h,m)::Constant
  V::Parameter
  Ψ::Form0

  ∂ₜ(Ψ) == ((-1 * (h^2)/(2*m))*Δ(Ψ) + V * Ψ) / (i*h)
end

Jordan_Kinderlehrer_Otto = @decapode begin
  (ρ,Ψ)::Form0
  β⁻¹::Constant
  ∂ₜ(ρ) == ∘(⋆,d,⋆)(d(Ψ)∧ρ) + β⁻¹*Δ(ρ)
end

Klausmeier_Eq2a = @decapode begin
  (n,w)::DualForm0
  dX::Form1
  (a,ν)::Constant
  ∂ₜ(w) == a - w - w * n^2 + ν * ℒ(dX, w)
end

Kealy = @decapode begin
  (n,w)::DualForm0
  dX::Form1
  (a,ν)::Constant
  ∂ₜ(w) == a - w - w * n^2 + ν * Δ(w)
end

Klausmeier_Eq2b = @decapode begin
  (n,w)::DualForm0
  m::Constant
  ∂ₜ(n) == w * n^2 - m*n + Δ(n)
end

Lejeune = @decapode begin
  ρ::Form0
  (μ, Λ, L)::Constant

  ∂ₜ(ρ) == ρ * (1 - μ + (Λ - 1) * ρ - ρ*ρ) +
    0.5 * (L*L - ρ) * Δ(ρ) - 0.125 * ρ * Δ(ρ) * Δ(ρ)
end

Turing_Continuous_Ring = @decapode begin
  (X, Y)::Form0
  (μ, ν, a, b, c, d)::Constant

  ∂ₜ(X) == a*X + b*Y + μ*Δ(X)
  ∂ₜ(Y) == c*X + d*Y + ν*Δ(X)
end

Mohamed_Equation10ForN2 = quote
  𝐮::Form1
  (P, 𝑝ᵈ)::Form0
  (negone, half, μ)::Constant

  ∂ₜ(𝐮) == 𝐮̇

  𝑝ᵈ == P + half * i(𝐮,𝐮)

  𝐮̇ == μ * ∘(d, ⋆, d, ⋆)(𝐮) + (negone)*⋆₁⁻¹(∧₁₀ₚᵈ(𝐮, ⋆(d(𝐮)))) + d(𝑝ᵈ)
end

end # Canon

