module Diagrams

using CombinatorialSpaces.ExteriorCalculus
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
import Catlab.Graphics: to_graphviz
using Catlab.Graphs

export eq_to_diagrams, FreeElDiagram, to_graphviz

@present ThElDiagram <: ThElements begin
  namee::Attr(El, Name)
  namea::Attr(Arr, Name)
end
@acset_type FreeElDiagram(ThElDiagram, index=[:src,:tgt, :dom, :cod])


function add_el!(ce, otype)
  o = vcat(incident(ce, [otype], :nameo)...)
  if isempty(o)
    o = [add_part!(ce, :Ob, nameo = otype)]
  end
  add_part!(ce, :El, namee = Symbol("anon_",gensym()), πₑ=first(o))
end

function add_els!(ce, obs)
  map(o -> add_el!(ce, o), obs)
end

function add_arr!(ce, htype, src, tgt)
  h = vcat(incident(ce, [htype], :nameh)...)
  if isempty(h)
    h = [add_part!(ce, :Hom, nameh = htype, dom=ce[src, :πₑ], cod=ce[tgt, :πₑ])]
  end
  add_part!(ce, :Arr, namea = Symbol(htype), src=src, tgt=tgt, πₐ = first(h))
end

function add_arrs!(ce, homs, srcs, tgts)
  map(i -> add_arr!(ce, homs[i], srcs[i], tgts[i]), 1:length(homs))
end

sp_otimes(expr, syntax) = expr isa syntax.Ob{:otimes} ? expr.args : [expr]

function m_proj_n(syntax, n, expr)
  #left = isempty(1:(n-1)) ? [] : syntax.proj2(foldl(⊗, expr.args[1:n-1]), foldl(⊗, expr.args[n:end]))
  #right = isempty((n+1):length(expr.args)) ? [] : syntax.proj1(expr.args[n], foldl(⊗, expr.args[(n+1):end]))
  #foldl(⋅, vcat(left, right))
  Symbol("π$n")
end

function add_term!(ce, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  if h == :compose
    cur_obs = in_els
    for arg in a
      cur_obs = add_term!(ce, cur_obs, arg, syntax)
    end
    cur_obs
  elseif h == :otimes
    new_obs = [add_els!(ce, ce[ce[in_els, :πₑ], :nameo]) for i in 1:length(a)]
    map(i -> add_arrs!(ce, syntax.mcopy.(ce[ce[in_els, :πₑ], :nameo]), in_els, new_obs[i]), 1:length(a))
    up_obs = vcat(map(i -> add_term!(ce, collect(new_obs[i]), a[i], syntax), 1:length(a))...)
  elseif h ∈ [:plus, :minus, :mult]
    new_obs = add_term!(ce, in_els, term, syntax; h=:otimes)
    op_types = sp_otimes(codom(first(a)), syntax)
    doub_types = op_types .⊗ op_types
    comb_els = add_els!(ce, doub_types)
    n_ops = length(op_types)
    map(n -> add_arrs!(ce, [Symbol("a", i + (n-1)*n_ops) for i in 1:n_ops], new_obs[(1:n_ops) .+ (n-1)*n_ops], comb_els), 1:length(a))
    out_obs = add_els!(ce, op_types)
    add_arrs!(ce, [syntax.plus(o) for o in op_types], comb_els, out_obs)
    out_obs
  elseif term isa syntax.Hom
    d_args = sp_otimes(dom(term), syntax)
    c_args = sp_otimes(codom(term), syntax)
    inp = first(in_els)
    if length(d_args) > 1
      inp = add_el!(ce, dom(term))
      add_arrs!(ce, [Symbol("a", i) for i in 1:length(in_els)], in_els, fill(inp, length(in_els)))
    end
    otp = add_el!(ce, codom(term))
    add_arr!(ce, term, inp, otp)
    otps = [otp]
    if length(c_args) > 1
      otps = add_els!(ce, c_args)
      add_arrs!(ce, [m_proj_n(syntax, i, ce[ce[otp, :πₑ], :nameo]) for i in 1:length(in_els)], fill(otp, length(otps)), otps)
    end
    otps
  else
    error("$h is not an implemented operator.")
  end
end

function rem_units!(ce, syntax)
  munits = findall(x -> x isa syntax.Ob{:munit}, ce[:nameo])
  unit_els = vcat(incident(ce, munits,:πₑ)...)
  if isempty(unit_els)
    return
  end
  i_wires = vcat(incident(ce, unit_els, :src)...)

  o_wires = vcat(incident(ce, unit_els, :tgt)...)

  ce[ce[i_wires, :tgt], :namee] .= ce[i_wires, :namea]
  ce[ce[o_wires, :src], :namee] .= ce[o_wires, :namea]

  rem_parts!(ce, :El, sort(unit_els))
  rem_parts!(ce, :Arr, sort(unique(vcat(i_wires, o_wires))))
end

function merge_nodes!(ce, args)
  source_node = first(args)

  i_wires = vcat(incident(ce, args, :src)...)
  o_wires = vcat(incident(ce, args, :tgt)...)

  set_subparts!(ce, i_wires, src = source_node)
  set_subparts!(ce, o_wires, tgt = source_node)

  rem_parts!(ce, :El, sort(unique(args[2:end])))
end

function merge_pairs!(ce, left, right)
  i_wires = incident(ce, right, :src)
  o_wires = incident(ce, right, :tgt)
  for i in 1:length(left)
    set_subparts!(ce, i_wires[i], src = left[i])
    set_subparts!(ce, o_wires[i], tgt = left[i])
  end
  rem_parts!(ce, :El, sort(unique(right)))
end

function merge_common!(ce)
  names = ce[:namee]
  for n in names
    vals = incident(ce, n, :namee)
    if length(vals) > 1
      merge_nodes!(ce, vals)
    end
  end
end
function eq_to_diagrams(d::Presentation)
  ce = FreeElDiagram{Any}()
  for eq in equations(d)
    first_ob = add_el!(ce, dom(eq[1]))
    left = add_term!(ce, [first_ob], eq[1], d.syntax)
    right = add_term!(ce, [first_ob], eq[2], d.syntax)
    merge_pairs!(ce, left, right)
  end
  rem_units!(ce, d.syntax)
  merge_common!(ce)
  ce
end

function to_graphviz(fd::FreeElDiagram; edge_len=1.9)
  nl(v) = (length(string(fd[v, :namee])) > 4 && string(fd[v, :namee])[1:4] == "anon") ? Dict(:shape => "point") : Dict(:label => "$(fd[v, :namee])")
	el(e) = Dict(:label => "$(fd[e, :namea])", :len => "$edge_len")
	g = Graph()
	migrate!(g, fd, Dict(:V=>:El, :E=>:Arr), Dict(:src=>:src, :tgt => :tgt))
	to_graphviz(PropertyGraph{Any}(g, nl, el; prog="neato"))
end

end
