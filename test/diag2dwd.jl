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

import Unicode
normalize_unicode(s::String) = Unicode.normalize(s, compose=true, stable=true, chartransform=Unicode.julia_chartransform)
normalize_unicode(s::Symbol)  = Symbol(normalize_unicode(String(s)))
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
# TODO: Add support for higher order functions.
#   - This is straightforward from a language perspective but unclear the best
#   - way to represent this in a Decapode ACSet.
@data Term begin
  Var(Symbol)
  Judge(Var, Symbol, Symbol) # Symbol 1: Form0 Symbol 2: X
  AppCirc1(Vector{Symbol}, Term)
  AppCirc2(Vector{Symbol}, Term, Term)
  App1(Symbol, Term)
  App2(Symbol, Term, Term)
  Plus(Term, Term)
  Tan(Term)
end

@data Equation begin
  Eq(Term, Term)
end

# A struct to store a complete Decapode
# TODO: Have the decopode macro compile to a DecaExpr
struct DecaExpr
  judgements::Vector{Judge}
  equations::Vector{Equation}
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
reduce_term!(t::Term, d::AbstractDecapode, syms::Dict{Symbol, Int}) =
  let ! = reduce_term!
    @match t begin
      Var(x) => syms[x]
      App1(f, t) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op1, src=!(t,d,syms), tgt=res_var, op1=f)
        return res_var
      end
      App2(f, t1, t2) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=!(t1,d,syms), proj2=!(t2,d,syms), res=res_var, op2=f)
        return res_var
      end
      AppCirc1(fs, t) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op1, src=!(t,d,syms), tgt=res_var, op1=fs)
        return res_var
      end
      AppCirc2(f, t1, t2) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=!(t1,d,syms), proj2=!(t2,d,syms), res=res_var, op2=fs)
        return res_var
      end
      Plus(t1, t2) => begin # TODO: plus is an Op2 so just fold it into App2
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=!(t1,d,syms), proj2=!(t2,d,syms), res=res_var, op2=:plus)
        return res_var
      end
      Tan(t) => begin 
        # TODO: this is creating a spurious variablbe with the same name
        txv = add_part!(d, :Var, type=:infer)
        tx = add_part!(d, :TVar, incl=txv)
        tanop = add_part!(d, :Op1, src=!(t,d,syms), tgt=txv, op1=DerivOp)
        return txv #syms[x._1]
      end
      _ => throw("Inline type judgements not yet supported!")
    end
  end

function eval_eq!(eq::Equation, d::AbstractDecapode, syms::Dict{Symbol, Int}) 
  @match eq begin
    Eq(t1, t2) => begin
      lhs_ref = reduce_term!(t1,d,syms)
      rhs_ref = reduce_term!(t2,d,syms)
      deletions = []
      # Make rhs_ref equal to lhs_ref and adjust all its incidents
      # Case rhs_ref is a Op1
      for rhs in incident(d, rhs_ref, :tgt)
        d[rhs, :tgt] = lhs_ref
        push!(deletions, rhs_ref)
      end
      # Case rhs_ref is a Op2
      for rhs in incident(d, rhs_ref, :res)
        d[rhs, :res] = lhs_ref
        push!(deletions, rhs_ref)
      end
      # TODO: delete unused vars. The only thing stopping me from doing 
      # this is I don't know if CSet deletion preserves incident relations
      rem_parts!(d, :Var, sort(deletions))
    end
  end
  return d
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
    eval_eq!(eq, d, symbol_table)
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
      eval_eq!(eq, d, symbol_table)
    end
    fill_names!(d)
    d[:name] = map(normalize_unicode,d[:name])
    return d
end

append_dot(s::Symbol) = Symbol(string(s)*'\U0307')

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
        d[t, :name] = append_dot(d[s,:name])
    end
    return d
end

function expand_operators(d::NamedDecapode)
  e = NamedDecapode{Symbol, Symbol, Symbol}()
  copy_parts!(e, d, (:Var, :TVar, :Op2))
  newvar = 0
  for op in parts(d, :Op1)
    if !isa(d[op,:op1], AbstractArray)
      add_part!(e, :Op1, op1=d[op,:op1], src=d[op, :src], tgt=d[op,:tgt])
    else
      for (i, step) in enumerate(d[op, :op1])
        if i == 1
          newvar = add_part!(e, :Var, type=:infer, name=Symbol("•_$(op)_$(i)"))
          add_part!(e, :Op1, op1=step, src=d[op, :src], tgt=newvar)
        elseif i == length(d[op, :op1])
          add_part!(e, :Op1, op1=step, src=newvar, tgt=d[op,:tgt])
        else
          newvar′ = add_part!(e, :Var, type=:infer, name=Symbol("•_$(op)_$(i)"))
          add_part!(e, :Op1, op1=step, src=newvar, tgt=newvar′)
          newvar = newvar′
        end
      end
    end
  end
  return e
end

abstract type AbstractCall end

struct UnaryCall <: AbstractCall 
    operator
    input
    output
end

Base.Expr(c::UnaryCall) = begin
    operator = c.operator
    if isa(c.operator, AbstractArray)
        operator = :(compose($operator))
    end
    return :($(c.output) = $operator($(c.input)))
end
                

struct BinaryCall <: AbstractCall 
    operator
    input1
    input2
    output
end

Base.Expr(c::BinaryCall) = begin
    operator = c.operator
    if isa(c.operator, AbstractArray)
        operator = :(compose($(c.operator)))
    end
    return :($(c.output) = $operator($(c.input1), $(c.input2)))
end

function get_vars_code(d::NamedDecapode, vars::Vector{Symbol})
    stmts = map(vars) do s
        ssymbl = QuoteNode(s)
        :($s = findnode(u, $ssymbl).values)
    end
    return quote $(stmts...) end
end


function set_tanvars_code(d::NamedDecapode, statevars::Vector{Symbol})
    tanvars = [(d[e, [:src,:name]], d[e, [:tgt,:name]]) for e in incident(d, :∂ₜ, :op1)]
    stmts = map(tanvars) do (s,t)
        ssymb = QuoteNode(s)
        :(findnode(du, $ssymb).values .= $t)
    end
    return quote $(stmts...) end
end


function compile(d::NamedDecapode, inputs::Vector)
    input_numbers = incident(d, inputs, :name)
    visited = falses(nparts(d, :Var))
    visited[collect(flatten(input_numbers))] .= true
    consumed1 = falses(nparts(d, :Op1))
    consumed2 = falses(nparts(d, :Op2))
    # FIXME: this is a quadratic implementation of topological_sort inlined in here.
    op_order = []
    for iter in 1:(nparts(d, :Op1) + nparts(d,:Op2))
        for op in parts(d, :Op1)
            s = d[op, :src]
            if !consumed1[op] && visited[s]
                # skip the derivative edges
                operator = d[op, :op1]
                t = d[op, :tgt]
                if operator == DerivOp
                    continue
                end
                consumed1[op] = true
                visited[t] = true
                sname = d[s, :name]
                tname = d[t, :name]
                c = UnaryCall(operator, sname, tname)
                push!(op_order, c)
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
                c = BinaryCall(operator, a1name, a2name, rname)
                push!(op_order, c)
            end
        end
    end
    assigns = map(Expr, op_order)
    ret = :(return)
    ret.args = d[d[:,:incl], :name]
    return quote f(du, u, p, t) = begin
        $(get_vars_code(d, inputs))
        $(assigns...)
        du .= 0.0
        $(set_tanvars_code(d, inputs))
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

# TODO: Write tests for recursive expressions

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

function compile_env(d::NamedDecapode)
  defs = quote end
  for op in d[:op1]
    if op == DerivOp
      continue
    end
    ops = QuoteNode(op)
    def = :($op = generate(mesh, $ops))
    push!(defs.args, def)
  end
  for op in d[:op2]
    if op == :+
      continue
    end
    ops = QuoteNode(op)
    def = :($op = generate(mesh, $ops))
    push!(defs.args, def)
  end
  return defs
end

function gensim(d::NamedDecapode, input_vars)
  d′ = expand_operators(d)
  defs = compile_env(d′)
  rhs = compile(d′, input_vars)
  quote
    function simulate(mesh)
      $defs
      return $rhs
    end
  end
end
## #DECAPODE Surface Syntax

# @data Term begin
#     Var(Symbol)
#     Judge(Var, Symbol, Symbol) # Symbol 1: Form0 Symbol 2: X
#     Eq(Term, Term)
#     AppCirc1(Vector{Symbol}, Var)
#     AppCirc2(Vector{Symbol}, Var, Var)
#     App1(Symbol, Var)
#     App2(Symbol, Var, Var)
#     Plus(Var, Var)
#     Tan(Var)
# end

term(s::Symbol) = Var(normalize_unicode(s))
term(expr::Expr) = begin
    @match expr begin
        Expr(a) => Var(normalize_unicode(a))
        Expr(:call, :∂ₜ, b) => Tan(term(b))
        Expr(:call, Expr(:call, :∘, a...), b) => AppCirc1(a, Var(b))
        Expr(:call, a, b) => App1(a, term(b))
        Expr(:call, Expr(:call, :∘, f...), x, y) => AppCirc1(f, Var(x), Var(y))
        Expr(:call, f, x, y) => App2(f, term(x), term(y))
        Expr(:call, :∘, a...) => (:AppCirc1, map(term, a))
        x => error("Cannot construct term from  $x")
    end
end

function parse_decapode(expr::Expr)
    stmts = map(expr.args) do line 
        @match line begin
            ::LineNumberNode => missing
            Expr(:(::), a, b) => Judge(Var(a),b.args[1], b.args[2])
            Expr(:call, :(==), lhs, rhs) => Eq(term(lhs), term(rhs))
            x => x
        end
    end |> skipmissing |> collect
    judges = []
    eqns = []
    for s in stmts
        if typeof(s) == Judge
            push!(judges, s)
        elseif typeof(s) == Eq
            push!(eqns, s)
        end
    end
    DecaExpr(judges, eqns)
end

term(:(∧₀₁(C,V)))

### Drawing of Decapodes

import Catlab.Graphics.Graphviz

#This could probably be made neater
function to_graphviz(d::NamedDecapode)::Graphviz.Graph
    #Similar to the to_graphviz in other implementations
    gv_name(v::Int) = "n$v"
    
    gv_path(e::Int) = [gv_name(d[:src][e]), gv_name(d[:tgt][e])]

    gp_name(p::Int) = "p$p"
    gp_proj1(p::Int) = [gp_name(p), gv_name(d[:proj1][p])]
    gp_proj2(p::Int) = [gp_name(p), gv_name(d[:proj2][p])]
    gp_projRes(p::Int) = [gp_name(p), gv_name(d[:res][p])]

    stmts = Graphviz.Statement[]

    reg_to_sub = Dict('0'=>'₀', '1'=>"₁", '2'=>'₂', '3'=>'₃', '4'=>'₄',
    '5'=>'₅', '6'=>'₆','7'=>'₇', '8'=>'₈', '9'=>'₉')

    toSub(digit::Char) = haskey(reg_to_sub, digit) ? reg_to_sub[digit] : digit

    #For variables, label grabs the stored variable name and its type and concatenate
    #label assumes dimension is single digit
    for v in parts(d, :Var)
        vertex_name = String(d[:name][v]) * ":Ω" * toSub(last(String(d[:type][v])))
        push!(stmts, Graphviz.Node(gv_name(v), Dict(:label=>vertex_name)))
    end

    #For unary ops, label mashes together all func symbol names into one string
    for e in parts(d, :Op1)
        #add composition symbols?
        edge_name = join(String.(d[:op1][e]))
        push!(stmts, Graphviz.Edge(gv_path(e), Dict(:label=>edge_name)))
    end

    #For binary ops, make temp product object, drop projections and drop result with op name
    for p in parts(d, :Op2)
        proj_space_name = "Ω" * toSub(last(String(d[:type][d[:proj1][p]]))) * "×" * "Ω" * toSub(last(String(d[:type][d[:proj2][p]])))
        push!(stmts, Graphviz.Node(gp_name(p), Dict(:label=>proj_space_name)))

        push!(stmts, Graphviz.Edge(gp_proj1(p), Dict(:label=>"proj₁", :style=>"dashed")))
        push!(stmts, Graphviz.Edge(gp_proj2(p), Dict(:label=>"proj₂", :style=>"dashed")))

        res_name = String(d[:op2][p])
        push!(stmts, Graphviz.Edge(gp_projRes(p), Dict(:label=>res_name)))
    end

    #Need to add user access for more customizability later
    Graphviz.Graph("G", true, "neato", stmts, Dict(), Dict(:shape=>"oval"), Dict())
end

using Test

@testset "Term Construction" begin
    @test term(:(Ċ)) == Var(:Ċ)
    @test_throws ErrorException term(:(∂ₜ{Form0}))
    @test term(Expr(:ϕ)) == Var(:ϕ)
    @test typeof(term(:(d₀(C)))) == App1
    @test typeof(term(:(∘(k, d₀)(C)))) == AppCirc1
    # @test term(:(∘(k, d₀)(C))) == AppCirc1([:k, :d₀], Var(:C)) #(:App1, ((:Circ, :k, :d₀), Var(:C)))
    # @test term(:(∘(k, d₀{X})(C))) == (:App1, ((:Circ, :k, :(d₀{X})), Var(:C)))
    @test_throws MethodError term(:(Ċ == ∘(⋆₀⁻¹{X}, dual_d₁{X}, ⋆₁{X})(ϕ)))
    @test term(:(∂ₜ(C))) == Tan(Var(:C))
    # @test term(:(∂ₜ{Form0}(C))) == App1(:Tan, Var(:C))
end

@testset "Recursive Expr" begin
  Recursion = quote
    x::Form0{X}
    y::Form0{X}
    z::Form0{X}

    ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
    y == F(f2(x), ρ(x,z))
  end

  recExpr = parse_decapode(Recursion)
  rdp = NamedDecapode(recExpr)
  show(rdp)

  @test nparts(rdp, :Var) == 9
  @test nparts(rdp, :TVar) == 1
  @test nparts(rdp, :Op1) == 5
  @test nparts(rdp, :Op2) == 3
end
Recursion = quote
  x::Form0{X}
  y::Form0{X}
  z::Form0{X}

  ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
  y == F(f2(x), ρ(x,z))
end

recExpr = parse_decapode(Recursion)
rdp = NamedDecapode(recExpr)

@testset "Diffusion Diagram" begin
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
    ddp = NamedDecapode(diffExpr)
    to_graphviz(ddp)

    @test nparts(ddp, :Var) == 3
    @test nparts(ddp, :TVar) == 1
    @test nparts(ddp, :Op1) == 3
    @test nparts(ddp, :Op2) == 0
end


@testset "Advection Diagram" begin
    Advection = quote
        C::Form0{X}
        V::Form1{X}
        ϕ::Form1{X}

        ϕ == ∧₀₁(C,V)
    end

    advdecexpr = parse_decapode(Advection)
    advdp = NamedDecapode(advdecexpr)
    @test nparts(advdp, :Var) == 3
    @test nparts(advdp, :TVar) == 0
    @test nparts(advdp, :Op1) == 0
    @test nparts(advdp, :Op2) == 1
end

@testset "Superposition Diagram" begin
    Superposition = quote
        C::Form0{X}
        Ċ::Form0{X}
        ϕ::Form1{X}
        ϕ₁::Form1{X}
        ϕ₂::Form1{X}

        ϕ == ϕ₁ + ϕ₂
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    superexp = parse_decapode(Superposition)
    supdp = NamedDecapode(superexp)
    @test nparts(supdp, :Var) == 5
    @test nparts(supdp, :TVar) == 1
    @test nparts(supdp, :Op1) == 2
    @test nparts(supdp, :Op2) == 1
end

@testset "AdvectionDiffusion Diagram" begin
    AdvDiff = quote
        C::Form0{X}
        Ċ::Form0{X}
        V::Form1{X}
        ϕ::Form1{X}
        ϕ₁::Form1{X}
        ϕ₂::Form1{X}
    
        # Fick's first law
        ϕ₁ ==  (k ∘ d₀)(C)
        # Diffusion equation
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ϕ₂ == ∧₀₁(C,V)

        ϕ == ϕ₁ + ϕ₂
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    advdiff = parse_decapode(AdvDiff)
    advdiffdp = NamedDecapode(advdiff)
    @test nparts(advdiffdp, :Var) == 6
    @test nparts(advdiffdp, :TVar) == 1
    @test nparts(advdiffdp, :Op1) == 4
    @test nparts(advdiffdp, :Op2) == 2
end


AdvDiff = quote
    C::Form0{X}
    Ċ::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}

    # Fick's first law
    ϕ₁ ==  (k ∘ d₀)(C)
    ϕ₂ == ∧₀₁(C,V)
    ϕ == ϕ₁ + ϕ₂
    # Diffusion equation
    ∂ₜ(C) == Ċ
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
end

advdiff = parse_decapode(AdvDiff)
advdiffdp = NamedDecapode(advdiff)
to_graphviz(advdiffdp)

compile(advdiffdp, [:C, :V])

