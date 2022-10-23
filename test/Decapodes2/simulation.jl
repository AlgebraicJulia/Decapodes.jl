using Decapodes
import Decapodes: compile, gensim, infer_states, infer_state_names

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics

using Test

@testset "Simulation Generation" begin

DiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ ==  ∘(k, d₀)(C)
    # Diffusion equation
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
add_scalar!(ddp, :k)
@test nparts(ddp, :Var) == 4

DiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}
    k::Scalar{ℝ}


    # Fick's first law
    ϕ ==  k * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
compile(expand_operators(ddp), [:C, :k])

@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k)
@test infer_state_names(ddp) == [:C, :k]

gensim(ddp)


DiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}
    k::Parameter{ℝ}


    # Fick's first law
    ϕ ==  k * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
compile(expand_operators(ddp), [:C, :k])


@test infer_state_names(ddp) == [:C, :k]
@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k(t))
gensim(ddp)
end
