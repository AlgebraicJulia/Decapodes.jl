### Drawing of Decapodes
import Catlab.Graphics.Graphviz
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs.PropertyGraphs

reg_to_sub = Dict('0'=>'₀', '1'=>"₁", '2'=>'₂', '3'=>'₃', '4'=>'₄',
    '5'=>'₅', '6'=>'₆','7'=>'₇', '8'=>'₈', '9'=>'₉', 'r'=>'•')

toSub(digit::Char) = get(reg_to_sub, digit, digit)

spacename(d, v) = begin
    t = d[v, :type]
    subscript = toSub(last(String(t)))
    dom = startswith(String(t), "Dual") ? "Ω̃" : "Ω"
    return "$dom$subscript"
end
varname(d, v) = "$(d[v, :name]):$(spacename(d, v))"

Graphics.to_graphviz(F::AbstractDecapode; kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(F; kw...))

decapode_edge_label(s::Symbol) = String(s)
decapode_edge_label(s::Vector{Symbol}) = join(String.(s), "⋅")


function Catlab.Graphics.to_graphviz_property_graph(d::AbstractNamedDecapode; kw...)
    pg = PropertyGraph{Any}(;kw...)
    vids = map(parts(d, :Var)) do v
      add_vertex!(pg, label=varname(d,v))
    end

    map(parts(d, :Op1)) do op
      s, t = d[op, :src], d[op, :tgt]
      add_edge!(pg, vids[s],vids[t], label=decapode_edge_label(d[op,:op1]))
    end

    map(parts(d, :Op2)) do op
      s, t = d[op, :proj1], d[op, :proj2]
      r = d[op, :res]
      v = add_vertex!(pg, label="$(spacename(d, s))×$(spacename(d,t))", shape="rectangle")
      add_edge!(pg, v, vids[s], label="π₁", style="dashed")
      add_edge!(pg, v, vids[t], label="π₂", style="dashed")
      add_edge!(pg, v, vids[r], label=decapode_edge_label(d[op, :op2]))
    end

    return pg
end

function Catlab.Graphics.to_graphviz_property_graph(d::SummationDecapode; kw...)
    tmp = NamedDecapode{Any, Any, Symbol}()
    # FIXME we have to copy to cast
    copy_parts!(tmp, d)
    G = to_graphviz_property_graph(tmp; kw...)
    findvid(G, d, v) = incident(G.graph, [Dict{Symbol, Any}(:label=>varname(d, v))], :vprops)
    white_nodes = map(parts(d, :Σ)) do s
        v = add_vertex!(G, label="Σ$s", shape="circle")
        u = d[s, :sum]
        matches = first(findvid(G, d, u))
        length(matches) == 1 || error("did not find a unique vertex match for Σ$s")
        uG = first(matches)
        add_edge!(G, v, uG, label="+")
        return v
    end
    for e in parts(d, :Summand)
        e = add_edge!(G, d[e, :summand], white_nodes[d[e, :summation]], style="dashed")
    end
    return G
end

savevizsvg(g, fname::String) = open(fname, "w") do fp
  run_graphviz(fp, to_graphviz(to_graphviz_property_graph(nsdp)), prog="neato", format="svg")
end