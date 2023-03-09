using Decapodes
import Decapodes: compile, gensim, infer_states, infer_state_names

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics

using Test

using MLStyle
using CombinatorialSpaces
using LinearAlgebra
using Distributions
using MultiScaleArrays

function generate(sd, my_symbol)
  op = @match my_symbol begin
    _ => default_dec_generate(sd, my_symbol, hodge)
    #= :⋆₀ => x->⋆(0,sd,hodge=DiagonalHodge())*x
    :⋆₁ => x->⋆(1, sd, hodge=DiagonalHodge())*x
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=DiagonalHodge())
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=DiagonalHodge())*x
    :d₀ => x->d(0,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y) =#
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length(args[1]))"); op(args...);end
  return (args...) ->  op(args...)
end

@testset "Simulation Generation" begin
DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ == ∘(k, d₀)(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
add_constant!(ddp, :k)
@test nparts(ddp, :Var) == 4

DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}
    k::Constant{ℝ}


    # Fick's first law
    ϕ == k * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
compile(expand_operators(ddp), [:C, :k])

@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k)
@test infer_state_names(ddp) == [:C, :k]

gensim(ddp)

torus = loadmesh(Torus_30x10())
c_dist = MvNormal([5, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in torus[:point]]

u₀ = construct(PhysicsState, [VectorForm(c)],Float64[], [:C])
du = construct(PhysicsState, [VectorForm(zero(c))],Float64[], [:C])

f = eval(gensim(expand_operators(ddp)))
fₘₛ = f(torus, generate)


DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}
    k::Parameter{ℝ}

    # Fick's first law
    ϕ == k * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
compile(expand_operators(ddp), [:C, :k])


@test infer_state_names(ddp) == [:C, :k]
@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k(t))
gensim(ddp)

f = eval(gensim(expand_operators(ddp)))
fₘₚ = f(torus, generate)

@test norm(fₘₛ(du, u₀, (k = _ -> 2.0,), 0)  - fₘₚ(du, u₀, (k=t->2.0,), 0)) < 1e-4

DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}
    k::Constant{ℝ}

    # Fick's first law
    ϕ == k * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
compile(ddp, [:C, :k])


gensim(ddp)

f = eval(gensim(ddp))
fₘₚ = f(torus, generate)

@test norm(fₘₛ(du, u₀, (k=2.0,), 0)  - fₘₚ(du, u₀, (k=2.0,), 0)) < 1e-4
 
# to solve the ODE over a duration, use the ODEProblem from OrdinaryDiffEq
# tₑ = 10
# using OrdinaryDiffEq
# prob = ODEProblem(fₘₛ,u₀, (0,tₑ), (k=2.0,))
# soln = solve(prob, Tsit5())
end
