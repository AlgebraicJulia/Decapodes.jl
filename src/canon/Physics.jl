module Physics

using DiagrammaticEquations
using DiagrammaticEquations.Deca
using ..Canon
using Markdown

@docapode("Mohamed Eq. 10, N2"
          ,""
          ,"A DEC formulation of fluid mechanics for surfaces. [mohamed_discrete_2016](@cite)"
          ,mohamed_flow
          ,begin
  (𝐮,w)::DualForm1
  (P, 𝑝ᵈ)::DualForm0
  μ::Constant

  𝑝ᵈ == P + 0.5 * ι₁₁(w,w)

  ∂ₜ(𝐮) == μ * ∘(d, ⋆, d, ⋆)(w) + (-1)*⋆₁⁻¹(∧ᵈᵖ₁₀(w, ⋆(d(w)))) + d(𝑝ᵈ)
end)

@docapode("Momentum"
  ,""
  ,""
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

@docapode("Ficks Law"
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

@docapode("Advection"
  ,"https://en.wikipedia.org/wiki/Advection"
  ,"Advection refers to the transport of a bulk along a vector field."
  ,advection,
  begin
    C::Form0
    (ϕ,V)::Form1

    ϕ == ∧₀₁(C,V)
  end
)

@docapode("Absorbed Shortwave Radiation"
  ,""
  ,"The proportion of light reflected by a surface is the **albedo**. The absorbed shortwave radiation is the complement of this quantity [budyko_effect_1969, sellers_global_1969](@cite)"
  ,absorbed_shortwave_radiation
  ,begin
    (Q, ASR)::Form0
    α::Constant

    ASR == (1 .- α) .* Q
  end
 )

@docapode("Outgoing Longwave Radiation"
  ,""
  ,"[budyko_effect_1969, sellers_global_1969](@cite)"
  ,:outgoing_longwave_radiation
  ,begin
    (Tₛ, OLR)::Form0
    (A,B)::Constant

    OLR == A .+ (B .* Tₛ)
  end
)

@docapode("Heat Transfer"
  ,""
  ,"[budyko_effect_1969](@cite) and [sellers_global_1969](@cite)"
  ,:heat_transfer
  ,begin
    (HT, Tₛ)::Form0
    (D,cosϕᵖ,cosϕᵈ)::Constant

    HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
  end
)

@docapode("Schroedinger"
  ,"https://en.wikipedia.org/wiki/Schrodinger_equation"
  ,"The evolution of the wave function over time."
  ,schroedinger
  ,begin
    (i,h,m)::Constant
    V::Parameter
    Ψ::Form0

    ∂ₜ(Ψ) == ((-1 * (h^2)/(2*m))*Δ(Ψ) + V * Ψ) / (i*h)
  end
)

@docapode("Navier-Stokes"
  ,"https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations"
  ,"Partial differential equations which describe the motion of viscous fluid surfaces."
  ,navier_stokes
  ,begin
    (V, V̇, G)::Form1{X}
    (ρ, ṗ, p)::Form0{X}

    V̇ == neg₁(L₁′(V, V)) +
        div₁(kᵥ(Δ₁(V) + third(d₀(δ₁(V)))), avg₀₁(ρ)) +
        d₀(half(i₁′(V, V))) +
        neg₁(div₁(d₀(p),avg₀₁(ρ))) +
        G
    ∂ₜ(V) == V̇
    ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
    ∂ₜ(p) == ṗ
  end
)

@docapode("Poiseuille"
  ,"https://en.wikipedia.org/wiki/Hagen-Poiseuille_equation"
  ,"A relation between the pressure drop in an incompressible and Newtownian fluid in laminar flow flowing through a long cylindrical pipe."
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

@docapode("Poiseuille Density"
  ,"https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation"
  ,""
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

@docapode("Jordan-Kinderlehrer-Otto"
  ,""
  ,"[jordan_variational_1998](@cite)"
  ,jko_scheme
  ,begin
    (ρ,Ψ)::Form0
    β⁻¹::Constant
    ∂ₜ(ρ) == ∘(⋆,d,⋆)(d(Ψ)∧ρ) + β⁻¹*Δ(ρ)
  end
)

@docapode("Oscillator"
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

@docapode(Lie
  ,"https://en.wikipedia.org/wiki/Lie_derivative"
  ,""
  ,lie
  ,begin
    C::Form0
    V::Form1
    dX::Form1

    V == ∘(⋆,⋆)(C ∧ dX)
  end
)

@docapode("Superposition"
  ,"https://en.wikipedia.org/wiki/Superposition"
  ,""
  ,superposition
  ,begin
    (C, Ċ)::Form0
    (ϕ, ϕ₁, ϕ₂)::Form1

    ϕ == ϕ₁ + ϕ₂
    Ċ == ⋆₀⁻¹(dual_d₁(⋆₁(ϕ)))
    ∂ₜ(C) == Ċ
  end
)

@docapode("IceBlockingWater"
  ,""
  ,"This model is an approach to modeling interface conditions of an obstruction within a flow developed by J. Fairbanks and L. Morris during the DARPA ASKEM hackathons. It should not be relied on as a realistic model of the physics of obstructions to flows."
  ,iceblockingwater
  ,begin
  h::Form0
  (𝐮,w)::DualForm1

  w == (1-σ(h)) ∧ᵖᵈ₀₁ 𝐮
  end
)
end
