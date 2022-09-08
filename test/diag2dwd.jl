using Decapodes
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Present
using Catlab.CategoricalAlgebra
# using Catlab.Graphics
# using Catlab.Syntax
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
# using CombinatorialSpaces.DiscreteExteriorCalculus: ∧
using LinearAlgebra

using Decapodes.Examples
using Decapodes.Diagrams
# using Decapodes.Simulations
# using Decapodes.Schedules

@present DiffusionSpace2D(FreeExtCalc2D) begin
  X::Space
  k::Hom(Form1(X), Form1(X)) # diffusivity of space, usually constant (scalar multiplication)
  proj₁_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
  proj₂_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
  sum₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
  prod₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
end


Diffusion = @decapode DiffusionSpace2D begin
    (C, Ċ₁, Ċ₂)::Form0{X}
    Ċ₁ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(k(d₀{X}(C)))))
    Ċ₂ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(d₀{X}(C))))
    ∂ₜ{Form0{X}}(C) == Ċ₁ + Ċ₂
end

# TODO: Replace the decapode macro with an ADT like this from MLStyle:
# Body => [Judge | Eqn]
# Judge => ([var])::Ob
# Eqn => Term == Term
# Term => Call(f, Term) | GENERATOR | Tan
# Tan => ∂ₜ(term)

# TODO: Write an interpreter for this ADT that generates databases over the following schema.

@present SchDecapode(FreeSchema) begin
    (Var, TVar, Op1, Op2)::Ob
    (Type, Operator)::AttrType
    # Name::AttrType
    src::Hom(Op1, Var)
    tgt::Hom(Op1, Var)
    proj1::Hom(Op2, Var)
    proj2::Hom(Op2, Var)
    res::Hom(Op2, Var)
    incl::Hom(TVar, Var)

    op1::Attr(Op1, Operator)
    op2::Attr(Op2, Operator)
    type::Attr(Var, Type)
    # name::Attr(Var, Name)
end


@abstract_acset_type AbstractDecapode 


@acset_type Decapode(SchDecapode,
  index=[:src, :tgt, :res, :incl, :op1, :op2, :type]) <: AbstractDecapode


begin
D = Diffusion
inputs = [D.ob_map[v] for v in vertices(J.graph)]
J = dom(D)
dp = Decapode{Any, Any}()

boxid = Vector{Int}()
for v in vertices(J.graph)
    vname = J.graph[v, :vname]
    Dv = D.ob_map[v]
    vid = add_part!(dp, :Var, type=(vname,Dv))
end
for e in edges(J.graph)
    De = D.hom_map[e]
    s = src(J.graph, e)
    t = tgt(J.graph, e)
    add_part!(dp, :Op1, src=s, tgt=t, op1=De)
    # push!(boxid, add_part!(dp, :Op1, src=Ds, tgt=Dt, operator=De))
end

end
dp
