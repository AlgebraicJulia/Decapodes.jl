module Diagrams

using CombinatorialSpaces.ExteriorCalculus
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
import Catlab.Graphics: to_graphviz
using Catlab.Graphs
using Catlab.Theories
using AlgebraicPetri
import AlgebraicPetri: Graph

export eq_to_diagrams, FreeElDiagram, to_graphviz, ElPetriNet, Graph

# Diagram of Elements from equations
####################################

@present ThElDiagram <: ThElements begin
  namee::Attr(El, Name)
  namea::Attr(Arr, Name)
end
@acset_type FreeElDiagram(ThElDiagram, index=[:src,:tgt, :dom, :cod])


function add_el!(ce::FreeElDiagram, otype)
  o = vcat(incident(ce, [otype], :nameo)...)
  if isempty(o)
    o = [add_part!(ce, :Ob, nameo = otype)]
  end
  add_part!(ce, :El, namee = Symbol("anon_",gensym()), πₑ=first(o))
end

function add_els!(ce::FreeElDiagram, obs)
  map(o -> add_el!(ce, o), obs)
end

function add_arr!(ce::FreeElDiagram, htype, src, tgt)
  h = vcat(incident(ce, [htype], :nameh)...)
  if isempty(h)
    h = [add_part!(ce, :Hom, nameh = htype, dom=ce[src, :πₑ], cod=ce[tgt, :πₑ])]
  end
  add_part!(ce, :Arr, namea = Symbol(htype), src=src, tgt=tgt, πₐ = first(h))
end

function add_arrs!(ce::FreeElDiagram, homs, srcs, tgts)
  map(i -> add_arr!(ce, homs[i], srcs[i], tgts[i]), 1:length(homs))
end

sp_otimes(expr, syntax) = expr isa syntax.Ob{:otimes} ? expr.args : [expr]

function m_proj_n(syntax, n, expr)
  #left = isempty(1:(n-1)) ? [] : syntax.proj2(foldl(⊗, expr.args[1:n-1]), foldl(⊗, expr.args[n:end]))
  #right = isempty((n+1):length(expr.args)) ? [] : syntax.proj1(expr.args[n], foldl(⊗, expr.args[(n+1):end]))
  #foldl(⋅, vcat(left, right))
  Symbol("π$n")
end

function add_term!(ce::FreeElDiagram, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  add_term!(Val{h}, ce, in_els, term, syntax; h=h, a=a)
end

function add_term!(::Type{Val{:compose}}, ce::FreeElDiagram, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  cur_obs = in_els
  for arg in a
    cur_obs = add_term!(ce, cur_obs, arg, syntax)
  end
  cur_obs
end

function add_term!(::Type{Val{:otimes}}, ce::FreeElDiagram, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
    new_obs = [add_els!(ce, ce[ce[in_els, :πₑ], :nameo]) for i in 1:length(a)]
    map(i -> add_arrs!(ce, syntax.mcopy.(ce[ce[in_els, :πₑ], :nameo]), in_els, new_obs[i]), 1:length(a))
    up_obs = vcat(map(i -> add_term!(ce, collect(new_obs[i]), a[i], syntax), 1:length(a))...)
end

function add_term!(::Type{Val{:plus}}, ce::FreeElDiagram, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
    new_obs = add_term!(ce, in_els, term, syntax; h=:otimes)
    op_types = sp_otimes(codom(first(a)), syntax)
    doub_types = op_types .⊗ op_types
    comb_els = add_els!(ce, doub_types)
    n_ops = length(op_types)
    map(n -> add_arrs!(ce, [Symbol("a", i + (n-1)*n_ops) for i in 1:n_ops], new_obs[(1:n_ops) .+ (n-1)*n_ops], comb_els), 1:length(a))
    out_obs = add_els!(ce, op_types)
    add_arrs!(ce, [syntax.plus(o) for o in op_types], comb_els, out_obs)
    out_obs
end

function add_term!(::Type{<:Val}, ce::FreeElDiagram, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
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
end



function rem_units!(ce::FreeElDiagram, syntax)
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

function merge_nodes!(ce::FreeElDiagram, args)
  source_node = first(args)

  i_wires = vcat(incident(ce, args, :src)...)
  o_wires = vcat(incident(ce, args, :tgt)...)

  set_subparts!(ce, i_wires, src = source_node)
  set_subparts!(ce, o_wires, tgt = source_node)

  rem_parts!(ce, :El, sort(unique(args[2:end])))
end

function merge_pairs!(ce::FreeElDiagram, left, right)
  i_wires = incident(ce, right, :src)
  o_wires = incident(ce, right, :tgt)
  for i in 1:length(left)
    set_subparts!(ce, i_wires[i], src = left[i])
    set_subparts!(ce, o_wires[i], tgt = left[i])
  end
  rem_parts!(ce, :El, sort(unique(right)))
end

function merge_common!(ce::FreeElDiagram)
  names = ce[:namee]
  for n in names
    vals = incident(ce, n, :namee)
    if length(vals) > 1
      merge_nodes!(ce, vals)
    end
  end
end

function to_graphviz(fd::FreeElDiagram; edge_len=1.9)
  nl(v) = (length(string(fd[v, :namee])) > 4 && string(fd[v, :namee])[1:4] == "anon") ? Dict(:shape => "point") : Dict(:label => "$(fd[v, :namee])")
	el(e) = Dict(:label => "$(fd[e, :namea])", :len => "$edge_len")
  g = Graphs.Graph()
	migrate!(g, fd, Dict(:V=>:El, :E=>:Arr), Dict(:src=>:src, :tgt => :tgt))
	to_graphviz(PropertyGraph{Any}(g, nl, el; prog="neato"))
end

# PetriNet of Elements
######################

# Need to add src/tgt information to Ob/Hom
@present ThElPetriNet(FreeSchema) begin
  Name::AttrType

  El::Ob
  Arr::Ob
  Src::Ob
  Tgt::Ob

  Ob::Ob
  Hom::Ob

  πₑ::Hom(El, Ob)
  πₐ::Hom(Arr, Hom)
  se::Hom(Src, El)
  sa::Hom(Src, Arr)
  te::Hom(Tgt, El)
  ta::Hom(Tgt, Arr)

  namee::Attr(El, Name)
  namea::Attr(Arr, Name)
  nameo::Attr(Ob, Name)
  nameh::Attr(Hom, Name)
end
@acset_type ElPetriNet(ThElPetriNet, index=[:se, :sa, :te, :ta])

function add_el!(ce::ElPetriNet, otype)
  o = vcat(incident(ce, [otype], :nameo)...)
  if isempty(o)
    o = [add_part!(ce, :Ob, nameo = otype)]
  end
  add_part!(ce, :El, namee = Symbol("anon_",gensym()), πₑ=first(o))
end

function add_els!(ce::ElPetriNet, obs)
  map(o -> add_el!(ce, o), obs)
end

function add_arr!(ce::ElPetriNet, htype, srcs, tgts)
  h = vcat(incident(ce, [htype], :nameh)...)
  if isempty(h)
    h = [add_part!(ce, :Hom, nameh = htype)]
  end
  arr = add_part!(ce, :Arr, namea = Symbol(htype), πₐ = first(h))
  add_parts!(ce, :Src, length(srcs), se=srcs, sa=fill(arr, length(srcs)))
  add_parts!(ce, :Tgt, length(tgts), te=tgts, ta=fill(arr, length(tgts)))
  arr
end

function add_arrs!(ce, homs, srcs, tgts)
  map(i -> add_arr!(ce, homs[i], srcs[i], tgts[i]), 1:length(homs))
end

function add_term!(ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  add_term!(Val{h}, ce, in_els, term, syntax; h=h, a=a)
end

function add_term!(::Type{Val{:compose}}, ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  cur_obs = in_els
  for arg in a
    cur_obs = add_term!(ce, cur_obs, arg, syntax)
  end
  cur_obs
end

function add_term!(::Type{Val{:otimes}}, ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  vcat(map(i -> add_term!(ce, in_els, a[i], syntax), 1:length(a))...)
end

function add_term!(::Type{Val{:oplus}}, ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
  up_obs = map(i -> add_term!(ce, in_els, a[i], syntax), 1:length(a))

  new_ob = add_el!(ce, foldl(⊕, [foldl(ExteriorCalculus.otimes, ce[ce[uo, :πₑ], :nameo]) for uo in up_obs]))
  add_arrs!(ce, [Symbol(:π, i) for i in 1:length(up_obs)], fill([new_ob], length(up_obs)), [u for u in up_obs])
  [new_ob]
end

function add_term!(::Type{Val{:plus}}, ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
    new_ob = add_term!(ce, in_els, term, syntax; h=:oplus)
    out_ob = add_el!(ce, codom(term))
    add_arr!(ce, syntax.plus(codom(term)), new_ob, [out_ob])
    [out_ob]
end

function add_term!(::Type{Val{:generator}}, ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
    d_args = sp_otimes(dom(term), syntax)
    c_args = sp_otimes(codom(term), syntax)

    codom_obs = add_els!(ce, c_args)
    add_arr!(ce, a[1], in_els, codom_obs)

    codom_obs
end

function add_term!(::Type{<:Val}, ce::ElPetriNet, in_els::Vector{Int}, term, syntax; h=head(term), a=args(term))
    d_args = sp_otimes(dom(term), syntax)
    c_args = sp_otimes(codom(term), syntax)

    codom_obs = add_els!(ce, c_args)
    add_arr!(ce, h, in_els, codom_obs)

    codom_obs
end



function rem_units!(ce::ElPetriNet, syntax)
  munits = findall(x -> x isa syntax.Ob{:munit}, ce[:nameo])
  unit_els = vcat(incident(ce, munits,:πₑ)...)
  if isempty(unit_els)
    return
  end
  i_wires = copy(ce[vcat(incident(ce, unit_els, :se)...), :sa])
  o_wires = copy(ce[vcat(incident(ce, unit_els, :te)...), :ta])

  i_els = ce[vcat(incident(ce, i_wires, :ta)...), :te]
  o_els = ce[vcat(incident(ce, o_wires, :sa)...), :se]

  ce[i_els, :namee] .= ce[i_wires, :namea]
  ce[o_els, :namee] .= ce[o_wires, :namea]

  rem_parts!(ce, :Src, sort(vcat(incident(ce, unit_els, :se)...)))
  rem_parts!(ce, :Tgt, sort(vcat(incident(ce, unit_els, :te)...)))
  rem_parts!(ce, :El, sort(unit_els))
  del_wires = unique(vcat(i_wires, o_wires))
  rem_parts!(ce, :Src, sort(vcat(incident(ce, del_wires, :sa)...)))
  rem_parts!(ce, :Tgt, sort(vcat(incident(ce, del_wires, :ta)...)))
  rem_parts!(ce, :Arr, sort(del_wires))
end

function merge_nodes!(ce::ElPetriNet, args)
  source_node = first(args)

  i_wires = vcat(incident(ce, args, :se)...)
  o_wires = vcat(incident(ce, args, :te)...)

  set_subparts!(ce, i_wires, se = source_node)
  set_subparts!(ce, o_wires, te = source_node)

  rem_parts!(ce, :El, sort(unique(args[2:end])))
end

function merge_pairs!(ce::ElPetriNet, left, right)

  i_wires = incident(ce, right, :se)
  o_wires = incident(ce, right, :te)
  for i in 1:length(left)
    set_subparts!(ce, i_wires[i], se = left[i])
    set_subparts!(ce, o_wires[i], te = left[i])
  end
  rem_parts!(ce, :El, sort(unique(right)))
end

function merge_common!(ce::ElPetriNet)
  names = ce[:namee]
  for n in names
    vals = incident(ce, n, :namee)
    if length(vals) > 1
      merge_nodes!(ce, vals)
    end
  end
end

function Graph(ep::ElPetriNet)
  pn = LabelledPetriNet()
  migrate!(pn, ep, Dict(:S => :El, :T => :Arr, :I => :Src, :O => :Tgt),
                   Dict(:is => :se, :os => :te, :it => :sa, :ot => :ta,
                        :tname => :namea, :sname => :namee))
  Graph(pn)
end

function eq_to_diagrams(d::Presentation; diagram=FreeElDiagram{Any}())
  ce = diagram
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

function to_monoidal(ep::ElPetriNet)
  ep2 = deepcopy(ep)

  proj_ind = findall(n -> first("$n") == 'π', ep[:namea])
  map(proj_ind) do pa
    proj_arrs = first(incident(ep, pa, :sa))
    el = ep[first(incident(ep, pa, :sa)), :se]
#    incident
  end
end

end
