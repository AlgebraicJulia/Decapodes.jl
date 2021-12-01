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

using ..Diagrams

export diag2dwd

sp_otimes(expr) = expr isa ObExpr{:otimes} ? expr.args : [expr]
n_args(expr) = length(sp_otimes(expr))

function eval_deps!(dwd, diagram, el, w2b, el2p; in_els = Dict{Int, Int}(), boundaries = [])
  if el in keys(el2p)
    return el2p[el]
  end

  inputs = incident(diagram, el, :tgt)

  labels = diagram[diagram[inputs, :πₐ], :nameh]
  el_types = sp_otimes(diagram[diagram[el, :πₑ], :nameo])
  outp = if length(inputs) > 1
    sort!(inputs, by = i->parse(Int, "$(diagram[diagram[i, :πₐ], :nameh])"[2:end]))
    vcat(map(w -> eval_deps!(dwd, diagram, diagram[w, :src], w2b, el2p; in_els = in_els), inputs)...)
  elseif length(inputs) == 0
    el ∈ keys(in_els) || error("Element $el has no dependencies, but is not defined in `in_els`")
    [(input_id(dwd), in_els[el])]
  else
    in_ps = eval_deps!(dwd, diagram, diagram[inputs[1], :src], w2b, el2p; in_els = in_els)
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
  bcs = incident(diagram, el, :src)
  filter!(a->"$(diagram[a, :namea])"[1] == '∂', bcs)
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

function diag2dwd(d; clean = false, calc_states = [])
  diagram = deepcopy(d)
  time_arrs = findall(h-> diagram[h, :nameh] isa HomExpr{:∂ₜ}, diagram[:πₐ])
  params = findall(e -> isempty(incident(diagram, e, :tgt)), 1:nparts(diagram, :El))

  state_vals = diagram[time_arrs, :src]
  in_els = unique(vcat(state_vals, params))
  out_els = copy(diagram[time_arrs, :tgt])
  out_names = [Symbol(:∂ₜ, diagram[e, :namee]) for e in state_vals]

  if !isempty(calc_states)
    calc_inds = findall(e->diagram[e, :namee] ∈ calc_states, state_vals)
    out_els = out_els[calc_inds]
    out_names = out_names[calc_inds]
  end

  in_types = [Dict(:name => diagram[e, :namee],
                   :type => diagram[diagram[e, :πₑ], :nameo])
              for e in in_els]
  out_types = [Dict(:name => out_names[e],
                    :type => diagram[diagram[out_els[e], :πₑ], :nameo])
               for e in 1:length(out_els)]
  rem_parts!(diagram, :Arr, time_arrs)

  dwd = WiringDiagram(collect(in_types), collect(out_types))
  w2b = map(1:nparts(diagram, :Arr)) do a
    w_type = diagram[diagram[a, :πₐ], :nameh]
    if w_type isa Symbol
      nothing
    elseif w_type isa HomExpr
      src_el = diagram[a, :src]
      tgt_el = diagram[a, :tgt]
      in_ports = sp_otimes(diagram[diagram[src_el, :πₑ], :nameo])
      out_ports = sp_otimes(diagram[diagram[tgt_el, :πₑ], :nameo])
      add_box!(dwd, Box(name(w_type), [Dict(:type=>ip) for ip in in_ports],
                        [Dict(:type=>op) for op in out_ports]))
    end
  end
  el2p = Dict{Int, Vector}()
  in_ps = vcat(map(out_els) do el
    eval_deps!(dwd, diagram, el, w2b, el2p; in_els=Dict(in_els[i] => i for i in 1:length(in_els)))
  end...)

  wires = map(enumerate(in_ps)) do (ip, op)
    Wire(out_types[ip], op, (output_id(dwd), ip))
  end
  add_wires!(dwd, wires)

  # Remove any boxes without any connecting wires
  if clean
    to_rem = Vector{Int64}()
    for b in 1:nparts(dwd.diagram, :Box)
      if isempty(vcat(incident(dwd.diagram, b, [:tgt, :in_port_box]),
                      incident(dwd.diagram, b, [:in_tgt, :in_port_box]))) &&
         isempty(vcat(incident(dwd.diagram, b, [:src, :out_port_box]),
                      incident(dwd.diagram, b, [:out_src, :out_port_box])))
         push!(to_rem, b)
      end
    end
    rem_boxes!(dwd, to_rem)
  end
  dwd
end
end
