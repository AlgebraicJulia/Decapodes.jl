""" Schedules

This module will contain the tooling which converts diagrams to DWDs. This will
also include optimization tooling relating to DWD schedules (compressing matrix
operations, pre-computing constant values, moving addition to reduce linear).
"""
module Schedules
using Catlab.Syntax
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Programs.DiagrammaticPrograms: NamedGraph

using ..Diagrams

export diag2dwd

sp_otimes(expr) = expr isa ObExpr{:otimes} ? expr.args : [expr]
n_args(expr) = length(sp_otimes(expr))

function eval_deps!(dwd, graph, el, w2b, el2p, obs; in_els = Dict{Int, Int}(), boundaries = [])
  if el in keys(el2p)
    return el2p[el]
  end

  inputs = incident(graph, el, :tgt)

  labels = graph[:ename]
  el_types = sp_otimes(obs[el])
  outp = if length(inputs) > 1
    sort!(inputs, by = i->parse(Int, "$(labels[i])"[2:end]))
    vcat(map(w -> eval_deps!(dwd, graph, graph[w, :src], w2b, el2p, obs; in_els = in_els), inputs)...)
  elseif length(inputs) == 0
    el ∈ keys(in_els) || error("Element $el has no dependencies, but is not defined in `in_els`")
    [(input_id(dwd), in_els[el])]
  else
    in_ps = eval_deps!(dwd, graph, graph[inputs[1], :src], w2b, el2p, obs; in_els = in_els)
    wires = map(enumerate(in_ps)) do (ip, op)
      Wire(port_value(dwd, Port(op[1], OutputPort, op[2])), op, (w2b[inputs[1]], ip))
    end
    add_wires!(dwd, wires)
    map(1:length(output_ports(dwd, w2b[inputs[1]]))) do p
      (w2b[inputs[1]], p)
    end
  end

  if el in keys(el2p)
    return el2p[el]
  end

  # Evaluate Boundary Conditions
  bcs = incident(graph, el, :src)
  filter!(a->"$(Symbol(graph[a, :ename]))"[1] == '∂', bcs)
  for bc in bcs
    wires = map(enumerate(outp)) do (ip, op)
      Wire(port_value(dwd, Port(op[1], OutputPort, op[2])), op, (w2b[bc], ip))
    end
    add_wires!(dwd, wires)
    outp .= map(1:length(output_ports(dwd, w2b[bc]))) do p
      (w2b[bc], p)
    end
  end
  el2p[el] = outp
  outp
end

name(a::HomExpr) = head(a) == :generator ? args(a)[1] : head(a)

function diag2dwd(diagram; clean = false, calc_states = [])
  homs = diagram.hom_map
  obs = diagram.ob_map
  graph = NamedGraph{Any, Any}()
  copy_parts!(graph, dom(diagram).graph)
  graph[:ename] .= homs

  # Expand homs which are composition
  for h in parts(graph, :E)
    h_name = graph[h, :ename]
    elsrc, eltgt = graph[h, :src], graph[h, :tgt]
    if h_name isa HomExpr{:compose}
      args = h_name.args
      rem_part!(graph, :E, h)
      verts = add_parts!(graph, :V, length(args) - 1, vname = :anon)
      append!(obs, codom.(args)[1:(end-1)])
      add_parts!(graph, :E, length(args), ename = args,
                          src = vcat([elsrc], verts),
                          tgt = vcat(verts, [eltgt]))
    end
  end

  pres = presentation(codom(diagram))
  time_arrs = findall(h-> h isa HomExpr{:∂ₜ}, graph[:ename])

  # TODO: Flip projection arrows for computation

  params = findall(e -> isempty(incident(graph, e, :tgt)), parts(graph, :V))

  state_vals = graph[time_arrs, :src]
  in_els = unique(vcat(state_vals, params))
  out_els = copy(graph[time_arrs, :tgt])
  out_names = [Symbol(:∂ₜ, graph[v, :vname]) for v in state_vals]

  if !isempty(calc_states)
    calc_inds = findall(v->graph[v, :vname] ∈ calc_states, state_vals)
    out_els = out_els[calc_inds]
    out_names = out_names[calc_inds]
  end

  in_types = [Dict(:name => graph[v, :vname],
                   :type => obs[v])
              for v in in_els]
  out_types = [Dict(:name => out_names[v],
                    :type => obs[v])
               for v in 1:length(out_els)]

  # FIXME: Hacky solution to ensuring time_arrs aren't included later
  # need more thoughtful approach to sorting through edges here
  rem_parts!(graph, :E, time_arrs)

  dwd = WiringDiagram(collect(in_types), collect(out_types))

  # Add all necessary boxes to the DWD based on arrows in the diagram
  w2b = map(parts(graph, :E)) do a
    w_type = graph[a, :ename]
    # Fix this if-case. This was meant for times when multiple arguments would
    # go into a single box
    if w_type isa Symbol
      nothing
    elseif w_type isa HomExpr
      src_el = graph[a, :src]
      tgt_el = graph[a, :tgt]
      in_ports = sp_otimes(obs[src_el])
      out_ports = sp_otimes(obs[tgt_el])
      add_box!(dwd, Box(name(w_type), [Dict(:type=>ip) for ip in in_ports],
                        [Dict(:type=>op) for op in out_ports]))
    end
  end

  # Mapping from elements (vertices in `graph`) to ports in the DWD
  # This ensures that elements which have already been computed are
  # not recomputed
  el2p = Dict{Int, Vector}()
  in_ps = vcat(map(out_els) do el
    eval_deps!(dwd, graph, el, w2b, el2p, obs; in_els=Dict(in_els[i] => i for i in 1:length(in_els)))
  end...)

  wires = map(enumerate(in_ps)) do (ip, op)
    Wire(out_types[ip], op, (output_id(dwd), ip))
  end
  add_wires!(dwd, wires)

  # Remove any boxes without any connecting wires
  if clean

  end

  dwd
end

end
