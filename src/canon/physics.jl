using DiagrammaticEquations
using DiagrammaticEquations.Deca

@docapode(Momentum
  ,"https://www.google.com"
  ,"desc"
  ,momentum
  ,begin
    (f,b)::Form0
    (v,V,g,Fᵥ,uˢ,v_up)::Form1
    τ::Form2
    U::Parameter

    uˢ̇ == ∂ₜ(uˢ)

    v_up == -1 * L(v,v) - L(V,v) - L(v,V) -
      f∧v - ∘(⋆,d,⋆)(uˢ)∧v - d(p) + b∧g - ∘(⋆,d,⋆)(τ) + uˢ̇ + Fᵥ

    uˢ̇ == force(U)
  
  end
)

@docapode(FicksLaw
  ,"https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion"
  ,"Equation for diffusion first stated by Adolf Fick. The diffusion flux is proportional to the concentration gradient."
  ,ficks_law
  ,begin
    C::Form0
    ϕ::Form1

    # Fick's first law
    ϕ ==  k(d₀(C))
  end
) 

@docapode(Advection
  ,"https://en.wikipedia.org/wiki/Advection"
  ,"Advection refers to the transport of a bulk along a vector field."
  ,advection,
  begin
    C::Form0
    (ϕ,V)::Form1

    ϕ == ∧₀₁(C,V)
  end
)

@docapode(AbsorbedShortwaveRadiation
  ,"https://www.google.com"
  ,"The proportion of light reflected by a surface is the **albedo**. The absorbed shortwave radiation is the complement of this quantity."
  ,absorbed_shortwave_radiation
  ,begin
    (Q, ASR)::Form0
    α::Constant

    ASR == (1 .- α) .* Q
  end
 )

# @docapode(:OutgoingLongwaveRadiation
#   ,"https://www.google.com"
#   ,"desc"
#   ,:outgoing_longwave_radiation
#   ,begin
#     (Tₛ, OLR)::Form0
#     (A,B)::Constant

#     OLR == A .+ (B .* Tₛ)
#   end
# )

# @docapode("Heat Transfer"
#   ,"https://www.google.com"
#   ,"desc"
#   ,"heat_transfer"
#   ,begin
#     (HT, Tₛ)::Form0
#     (D,cosϕᵖ,cosϕᵈ)::Constant

#     HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
#   end
# )

@docapode(Schoedinger
  ,"https://en.wikipedia.org/wiki/Schrodinger_equation"
  ,"The evolution of the wave functioin over time."
  ,schroedinger
  ,begin
    (i,h,m)::Constant
    V::Parameter
    Ψ::Form0

    ∂ₜ(Ψ) == ((-1 * (h^2)/(2*m))*Δ(Ψ) + V * Ψ) / (i*h)
  end
)

@docapode(NavierStokes
  ,"https://en.wikipedia.org/wiki/Navier_Stokes_equation"
  ,"Partial differential equations which describe the motion of viscous fluid surfaces."
  ,navier_stokes
  ,begin
    (V, V̇, G)::Form1{X}
    (ρ, ṗ, p)::Form0{X}

    # TODO: Find the right LHS for the next line
    V == Δ₁(V) + third(d₀(δ₁(V)))

    ∂ₜ(V) == neg₁(L₁′(V, V)) + 
      kᵥ(Δ(V) + (1/3)*d(δ(V))) / avg₀₁(ρ) +
      d(0.5 * (i₁′(V, V))) +
      -1 * d(p) / avg₀₁(ρ) +
      G

    ∂ₜ(p) == -1 * ⋆(L₀(V, ⋆(p)))
  end
),

@docapode(Poiseuille
  ,"https://www.google.com"
  ,"desc"
  ,poiseuille
  ,begin
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
)

@docapode(PoiseuilleDensity
  ,"https://en.wikipedia.org/wiki/hagen-poiseuille_density"
  ,"desc"
  ,poiseuille_density
  ,begin
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
)

@docapode(JordanKinderlehrerOtto
  ,"https://www.google.com"
  ,"desc"
  ,jko_scheme
  ,begin
    (ρ,Ψ)::Form0
    β⁻¹::Constant
    ∂ₜ(ρ) == ∘(⋆,d,⋆)(d(Ψ)∧ρ) + β⁻¹*Δ(ρ)
  end
)

@docapode(Oscillator
  ,"https://en.wikipedia.org/wiki/Harmonic_oscillator"
  ,"Equation governing the motion of an object whose acceleration is negatively-proportional to its position."
  ,oscillator
  ,begin
    X::Form0
    V::Form0
    k::Constant

    ∂ₜ(X) == V
    ∂ₜ(V) == -k*(X)
  end
)

# TODO
@docapode(Lie
  ,"https://en.wikipedia.org/wiki/lie_derivative"
  ,"..."
  ,lie
  ,begin
    C::Form0
    V::Form1
    dX::Form1

    V == ∘(⋆,⋆)(C ∧ dX)
  end
)

@docapode(Superposition
  ,"https://en.wikipedia.org/wiki/superposition"
  ,"desc"
  ,superposition
  ,begin
    (C, Ċ)::Form0
    (ϕ, ϕ₁, ϕ₂)::Form1

    ϕ == ϕ₁ + ϕ₂
    Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
    ∂ₜ(C) == Ċ
  end
)




















# Diffusion = @decapode DiffusionQuantities begin
#   (C, Ċ)::Form0{X}
#   ϕ::Form1{X}

#   # Fick's first law
#   ϕ ==  k(d₀{X}(C))
#   # Diffusion equation
#   Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
#   ∂ₜ{Form0{X}}(C) == Ċ
# end

# Diffusion = @decapode DiffusionQuantities begin
#   C::Form0{X}
#   ϕ::Form1{X}

#   # Fick's first law
#   ϕ ==  k(d₀{X}(C))
# end

# Advection = @decapode DiffusionQuantities begin
#   C::Form0{X}
#   (V, ϕ)::Form1{X}
#   ϕ == ∧₀₁{X}(C,V)
# end

# Superposition = @decapode DiffusionQuantities begin
#   (C, Ċ)::Form0{X}
#   (ϕ, ϕ₁, ϕ₂)::Form1{X}

#   ϕ == ϕ₁ + ϕ₂
#   Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
#   ∂ₜ{Form0{X}}(C) == Ċ
# end

