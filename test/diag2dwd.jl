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
using Base.Iterators

using Decapodes.Examples
using Decapodes.Diagrams
# using Decapodes.Simulations
# using Decapodes.Schedules

DerivOp = Symbol("∂ₜ")

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
    src::Hom(Op1, Var)
    tgt::Hom(Op1, Var)
    proj1::Hom(Op2, Var)
    proj2::Hom(Op2, Var)
    res::Hom(Op2, Var)
    incl::Hom(TVar, Var)
    
    op1::Attr(Op1, Operator)
    op2::Attr(Op2, Operator)
    type::Attr(Var, Type)
end

@present SchNamedDecapode <: SchDecapode begin
    Name::AttrType
    name::Attr(Var, Name)
end

@abstract_acset_type AbstractDecapode

@acset_type Decapode(SchDecapode,
  index=[:src, :tgt, :res, :incl, :op1, :op2, :type]) <: AbstractDecapode

@acset_type NamedDecapode(SchNamedDecapode,
  index=[:src, :tgt, :res, :incl, :op1, :op2, :type, :name]) <: AbstractDecapode

# to_decapode helper functions
reduce_lhs!(eq::Eq, d::AbstractDecapode, syms::Dict{Symbol, Int}) =
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
        txv = add_part!(d, :Var, type=:infer)
        tx = add_part!(d, :TVar, incl=txv)
        tanop = add_part!(d, :Op1, src=syms[x._1], tgt=txv, op1=DerivOp)
        return txv #syms[x._1]
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
reduce_rhs!(eq::Eq, d::AbstractDecapode, syms::Dict{Symbol, Int}, lhs_ref::Int) =
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
function Decapode(e::DecaExpr)
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

function NamedDecapode(e::DecaExpr)
    d = NamedDecapode{Any, Any, Symbol}()
    symbol_table = Dict{Symbol, Int}()
    for judgement in e.judgements
      var_id = add_part!(d, :Var, name=judgement._1._1, type=judgement._2)
      symbol_table[judgement._1._1] = var_id
    end
    for eq in e.equations
      v = reduce_lhs!(eq, d, symbol_table)
      reduce_rhs!(eq, d, symbol_table, v)
    end
    fill_names!(d)
    return d
end
"""    fill_names!

add new variable names to all the variables that don't have names.
"""
function fill_names!(d::NamedDecapode)
    bulletcount = 1
    for i in parts(d, :Var)
        if !isassigned(d[:,:name],i)
            d[i,:name] = Symbol("•$bulletcount")
            bulletcount += 1
        end
    end
    for e in incident(d, :∂ₜ, :op1)
        s = d[e,:src]
        t = d[e, :tgt]
        @show (e, s, t)
        d[t, :name] = Symbol("$(d[s,:name])̇")
    end
    return d
end

function compile(d::NamedDecapode, inputs::Vector)
    input_tuple = :(())
    append!(input_tuple.args, inputs)
    input_numbers = incident(d, inputs, :name)
    visited = falses(nparts(d, :Var))
    visited[collect(flatten(input_numbers))] .= true
    consumed1 = falses(nparts(d, :Op1))
    consumed2 = falses(nparts(d, :Op2))
    assigns = Expr[]
    # FIXME: this is a quadratic implementation of topological_sort inlined in here.
    for iter in 1:(nparts(d, :Op1) + nparts(d,:Op2))
        for op in parts(d, :Op1)
            s = d[op, :src]
            if !consumed1[op] && visited[s]
                t = d[op, :tgt]
                sname = d[s, :name]
                tname = d[t, :name]
                operator = d[op, :op1]
                # skip the derivative edges
                if operator == DerivOp
                    continue
                end
                consumed1[op] = true
                visited[t] = true
                if isa(operator, AbstractArray)
                    operator = :(compose($operator))
                end
                assignment = :($tname = $operator($sname))
                push!(assigns, assignment)
            end
        end

        for op in parts(d, :Op2)
            arg1 = d[op, :proj1]
            arg2 = d[op, :proj2]
            if !consumed2[op] && visited[arg1] && visited[arg2]
                r = d[op, :res]
                a1name = d[arg1, :name]
                a2name = d[arg2, :name]
                rname  = d[r, :name]
                operator = d[op, :op2]
                consumed2[op] = true
                visited[r] = true
                if isa(operator, AbstractArray)
                    operator = :(compose($operator))
                end
                assignment = :($rname = $operator($a1name, $a2name))
                push!(assigns, assignment)
            end
        end
    end
    ret = :(return)
    ret.args = d[d[:,:incl], :name]
    return quote f($(input_tuple)) = begin
        $(assigns...)
        $ret
    end; end
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
diffusion_cset = Decapode(diffusion_d)
diffusion_cset_named = NamedDecapode(diffusion_d)
# A test with expressions on LHS (i.e. temporary variables must be made)
# TODO: we should be intelligent and realize that the LHS of the first two
# equations are the same and so can share a new variable
eqs = [Eq(Plus(Var(:Ċ₁), Var(:Ċ₂)), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Plus(Var(:Ċ₁), Var(:Ċ₂)), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus(Var(:Ċ₁), Var(:Ċ₂)))    
]
test_d = DecaExpr(js, eqs)
test_cset = Decapode(test_d)
test_cset_named = NamedDecapode(test_d)

all(isassigned(test_cset_named[:name], i) for i in parts(test_cset_named,:Var))

sup_js = js = [Judge(Var(:C), :Form0, :X), 
Judge(Var(:ϕ₁), :Form0, :X),
Judge(Var(:ϕ₂), :Form0, :X)
]
sup_eqs = [Eq(Var(:ϕ₁), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Var(:ϕ₂), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus(Var(:ϕ₁), Var(:ϕ₂)))    
]
sup_d = DecaExpr(sup_js, sup_eqs)
sup_cset = Decapode(sup_d)
sup_cset_named = NamedDecapode(sup_d)


compile(diffusion_cset_named, [:C,])
compile(test_cset_named, [:C,])
compile(sup_cset_named, [:C,])

