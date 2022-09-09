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
using MLStyle

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

# Definition of the Decapodes DSL AST
# TODO: do functions and tans need to be parameterized by a space?
# TODO: make recursive
@data Term begin
  Var(Symbol)
  Judge(Var, Symbol, Symbol) # Symbol 1: Form0 Symbol 2: X
  Eq(Term, Term)
  AppCirc1(Vector{Symbol}, Var)
  AppCirc2(Vector{Symbol}, Var, Var)
  App1(Symbol, Var)
  App2(Symbol, Var, Var)
  Plus(Var, Var)
  Tan(Var)
end

# A struct to store a complete Decapode
# TODO: Have the decopode macro compile to a DecaExpr
struct DecaExpr
  judgements::Vector{Judge}
  equations::Vector{Eq}
end

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

# to_decapode helper functions
reduce_lhs!(eq::Eq, d::Decapode, syms::Dict{Symbol, Int}) =
  let ! = reduce_lhs! # This will be needed once we upgrade to a recursive grammar
    @match eq._1 begin
      Var(x) => syms[x]
      App1(f, x) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op1, src=syms[x._1], tgt=res_var, op1=f)
        return res_var
      end
      App2(f, x, y) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=res_var, op2=f)
        return res_var
      end
      AppCirc1(fs, x) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op1, src=syms[x._1], tgt=res_var, op1=fs)
        return res_var
      end
      AppCirc2(fs, x, y) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=res_var, op2=fs)
        return res_var
      end
      Tan(x) => begin
        tx = add_part!(d, :TVar, incl=syms[x._1])
        return syms[x._1]
      end
      Plus(x, y) => begin # TODO: plus is an Op2 so just fold it into App2
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=res_var, op2=:plus)
        return res_var
      end
      _ => -1 # TODO: make this throw an error or something
    end
  end
# TODO: lots of code duplication between reduce_lhs! and reduce_rhs!
# The duplicate code should be abstracted into another helper function
reduce_rhs!(eq::Eq, d::Decapode, syms::Dict{Symbol, Int}, lhs_ref::Int) =
  let ! = reduce_rhs! # Again only necessary once we upgrade language
    @match eq._2 begin
      App1(f, x) => add_part!(d, :Op1, src=syms[x._1], tgt=lhs_ref, op1=f)
      App2(f, x, y) => add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=lhs_ref, op2=f)
      AppCirc1(fs, x) => add_part!(d, :Op1, src=syms[x._1], tgt=lhs_ref, op1=fs)
      AppCirc2(fs, x, y) => add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=lhs_ref, op2=fs)
      Plus(x, y) => add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=lhs_ref, op2=:plus)
      _ => -1 # TODO: Throw an error or handle case where RHS is a raw variable or tangent
    end
  end

""" Takes a DecaExpr (i.e. what should be constructed using the @decapode macro)
and gives a Decapode ACSet which represents equalities as two operations with the
same tgt or res map.
"""
# Just to get up and running, I tagged untyped variables with :infer
# TODO: write a type checking/inference step for this function to 
# overwrite the :infer tags
function to_decapode(e::DecaExpr)::Decapode
  d = Decapode{Any, Any}()
  symbol_table = Dict{Symbol, Int}()
  for judgement in e.judgements
    var_id = add_part!(d, :Var, type=(judgement._2, judgement._3))
    symbol_table[judgement._1._1] = var_id
  end
  for eq in e.equations
    v = reduce_lhs!(eq, d, symbol_table)
    reduce_rhs!(eq, d, symbol_table, v)
  end
  return d
end

# Tests
#######

# Construct roughly what the @decapode macro should return for Diffusion
js = [Judge(Var(:C), :Form0, :X), 
      Judge(Var(:Ċ₁), :Form0, :X),
      Judge(Var(:Ċ₂), :Form0, :X)
]
# TODO: Do we need to handle the fact that all the functions are parameterized by a space?
eqs = [Eq(Var(:Ċ₁), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Var(:Ċ₂), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus(Var(:Ċ₁), Var(:Ċ₂)))
]
diffusion_d = DecaExpr(js, eqs)
diffusion_cset = to_decapode(diffusion_d)

# A test with expressions on LHS (i.e. temporary variables must be made)
# TODO: we should be intelligent and realize that the LHS of the first two
# equations are the same and so can share a new variable
eqs = [Eq(Plus(Var(:Ċ₁), Var(:Ċ₂)), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Plus(Var(:Ċ₁), Var(:Ċ₂)), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus(Var(:Ċ₁), Var(:Ċ₂)))    
]
test_d = DecaExpr(js, eqs)
test_cset = to_decapode(test_d)

#=begin
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
dp=#
