""" Generating Decapodes from Petri nets

This module provides support for generating a simulatable Decapode from Petri
nets generated in the AlgebraicPetri library.
"""

module PetriNets

using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphics, Catlab.Programs
using CombinatorialSpaces.ExteriorCalculus
using Catlab.Programs.DiagrammaticPrograms: NamedGraph
using ..Simulations
using ..AlgebraicPetri
using LinearAlgebra
using Unicode

export pn2dec, expand_pres!, gen_functions

i2sub = Dict(
  '0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', '5'=>'₅',
  '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉', '-'=>'₋'
)
i2sup = Dict(
  '0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', '5'=>'⁵',
  '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻'
)

function neg_term(graph, obs, homs, prs, val::Int64)
  res = add_part!(graph, :V, vname=Symbol("$(gensym())"))
  push!(obs, prs.syntax.Form0(prs[:X]))
  add_part!(graph, :E, ename=nothing, src=val, tgt=res)
  push!(homs, prs[Symbol("neg₀")])
  return res
end

function scale_term(graph, obs, homs, prs, val::Int64, sc::Int64)
  res = add_part!(graph, :V, vname="$(gensym())")
  push!(obs, prs.syntax.Form0(prs[:X]))
  add_part!(graph, :E, ename=nothing, src=val, tgt=res)
  push!(homs, prs[Symbol("mult_$(sc)")])
  return res
end

function exp_term(graph, obs, homs, prs, val::Int64, sc::Int64)
  res = add_part!(graph, :V, vname=Symbol("$(gensym())"))
  push!(obs, prs.syntax.Form0(prs[:X]))
  add_part!(graph, :E, ename=nothing, src=val, tgt=res)
  push!(homs, prs[Symbol("^$(sc)")])
  return res
end

function rate_term(graph, obs, homs, prs, val::Int64, tr_name::Symbol, tgt::Int64)
  add_part!(graph, :E, ename=nothing, src=val, tgt=tgt)
  push!(homs, prs[Symbol("k_", tr_name)])
  return tgt
end

function rate_term(graph, obs, homs, prs, val::Int64, tr_name::Symbol)
  res = add_part!(graph, :V, vname=Symbol("$(gensym())"))
  push!(obs, prs.syntax.Form0(prs[:X]))
  rate_term(graph, obs, homs, prs, val, tr_name, res)
end

function sum_dec!(graph, obs, homs, prs, vals::Int64...)
  bin_exp!(graph, obs, homs, prs, :sum₀, vals...)
end

function prod_dec!(graph, obs, homs, prs, vals::Int64...)
  bin_exp!(graph, obs, homs, prs, :prod₀, vals...)
end

function bin_exp!(graph, obs, homs, prs, exp, vals::Int64...)
  if length(vals) == 1
    first(vals)
  else
    form0 = prs.syntax.Form0(prs[:X])
    form0x0 = prs.syntax.otimes(form0, form0)
    cur_val = first(vals)
    bin_int = add_parts!(graph, :V, length(vals)-1,
              vname=[Symbol("$(gensym())") for i in 1:(length(vals)-1)])
    bin_res = add_parts!(graph, :V, length(vals)-1,
              vname=[Symbol("$(gensym())") for i in 1:(length(vals)-1)])
    append!(obs, fill(form0x0, length(vals)-1))
    append!(obs, fill(form0, length(vals)-1))
    for i in 1:(length(vals)-1)
      add_parts!(graph, :E, 3,
                ename =[nothing, nothing, nothing],
                tgt=[vals[i+1], cur_val, bin_res[i]],
                src=[bin_int[i], bin_int[i], bin_int[i]])
      append!(homs, [prs[:proj₁_⁰⁰₀], prs[:proj₂_⁰⁰₀], prs[exp]])
      cur_val = bin_res[i]
    end
    cur_val
  end
end

""" pn2dec(prs::Presentation, pn::LabelledReactionNet)

Generates a Decapode diagram which represents the law of mass action applied to
the petri net `pn` with the syntax from the presentation `prs`.
"""
function pn2dec(prs, pn::LabelledReactionNet)
  synt = prs.syntax

  # Initial state variables
  graph = NamedGraph{Symbol,Union{Symbol,Nothing}}()
  svars = add_parts!(graph, :V, ns(pn), vname=snames(pn))
  tvars = add_parts!(graph, :V, nt(pn), vname=tnames(pn))
  obs  = Vector{Any}(fill(synt.Form0(prs[:X]), ns(pn) + nt(pn))) #"($(join(vcat(snames(pn), tnames(pn)), ", ")))::Form0{X}"
  homs = Vector{Any}()
  map(1:nt(pn)) do t
    inputs = pn[incident(pn, t, :it), :is]
    res = prod_dec!(graph, obs, homs, prs, svars[inputs]...)
    rate_term(graph, obs, homs, prs, res, tname(pn, t), tvars[t])
  end

  res_s = map(1:ns(pn)) do s
    inps = tvars[pn[incident(pn, s, :os), :ot]]
    otps = tvars[pn[incident(pn, s, :is), :it]]
    otp_scaled = [neg_term(graph, obs, homs, prs, o) for o in otps]
    if isempty(vcat(inps, otp_scaled))
      res = add_part!(graph, :V, vname=Symbol("$(gensym())"))
      push!(obs, prs.syntax.Form0(prs[:X]))
      res
    else
      sum_dec!(graph, obs, homs, prs, vcat(inps, otp_scaled)...)
    end
  end

  set_subparts!(graph, res_s,
    vname=[length("$name") > 1 ?
            Symbol("∂ₜ", "$name") : Symbol(Unicode.normalize("$name"* '̇'))
            for name in snames(pn)])
  FinFunctor(obs, homs, FinCat(graph), FinCat(prs))
end

""" expand_pres!(pres::Presentation, pn::LabelledReactionNet)

Expands the syntax of the presentation `pres` to include the necessary operators for expressing mass action on the Petri net `pn`.
"""
function expand_pres!(pres, pn::LabelledReactionNet)
  synt = pres.syntax
  form0 = synt.Form0(pres[:X])
  form0x0 = synt.otimes(form0, form0)
  gens = Dict(:neg₀ => synt.Hom(:neg₀, form0, form0),
              :sum₀ => synt.Hom(:sum₀, form0x0, form0),
              :prod₀ => synt.Hom(:prod₀, form0x0, form0),
              :proj₁_⁰⁰₀ => synt.Hom(:proj₁_⁰⁰₀, form0x0, form0),
              :proj₂_⁰⁰₀ => synt.Hom(:proj₂_⁰⁰₀, form0x0, form0))
  for t in tnames(pn)
    gens[t] = synt.Hom(Symbol("k_$t"), form0, form0)
  end

  for (k,v) in gens
    if !has_generator(pres, k)
      add_generator!(pres, v)
    end
  end
end

""" gen_functions(pn::LabelledReactionNet, s::EmbeddedDeltaDualComplex2D)

Generates the functions necessary for executing the Petri net after it has been
converted to a Decapode.
"""
function gen_functions(pn::LabelledReactionNet, s)
  funcs = Dict{Symbol, Any}()
  funcs[:sum₀] = Dict(:operator => (s,a,b) -> (s .= a .+ b), :type => InPlaceFunc())
  funcs[:prod₀] = Dict(:operator => (s,a,b) -> (s .= a .* b), :type => InPlaceFunc())
  funcs[:neg₀] = Dict(:operator => (s,a) -> (s .= -1 * a), :type => InPlaceFunc())
  for t in 1:nt(pn)
    funcs[Symbol("k_", tname(pn, t))] = Dict(:operator => rate(pn, t) * I(nv(s)), :type => MatrixFunc())
  end
  funcs
end

end
