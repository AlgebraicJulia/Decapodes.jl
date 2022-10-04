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


#This could probably be made neater
function Catlab.Graphics.to_graphviz(d::AbstractNamedDecapode)::Graphviz.Graph
    #Similar to the to_graphviz in other implementations
    gv_name(v::Int) = "n$v"
    
    gv_path(e::Int) = [gv_name(d[e, :src]), gv_name(d[e, :tgt])]

    gp_name(p::Int) = "p$p"
    gp_proj1(p::Int) = [gp_name(p), gv_name(d[p, :proj1])]
    gp_proj2(p::Int) = [gp_name(p), gv_name(d[p, :proj2])]
    gp_projRes(p::Int) = [gp_name(p), gv_name(d[p, :res])]

    stmts = Graphviz.Statement[]

    #For variables, label grabs the stored variable name and its type and concatenate
    #label assumes dimension is single digit
    for v in parts(d, :Var)
        vertex_name = varname(d, v)
        push!(stmts, Graphviz.Node(gv_name(v), Dict(:label=>vertex_name)))
    end

    #For unary ops, label mashes together all func symbol names into one string
    for e in parts(d, :Op1)
        #add composition symbols?
        edge_name = join(String.(d[e, :op1]))
        push!(stmts, Graphviz.Edge(gv_path(e), Dict(:label=>edge_name)))
    end

    #For binary ops, make temp product object, drop projections and drop result with op name
    for p in parts(d, :Op2)
        v1 = d[p, :proj1]
        v2 = d[p, :proj2]
        proj_space_name = "$(spacename(d, v1)) × $(spacename(d, v2))"

        push!(stmts, Graphviz.Node(gp_name(p), Dict(:label=>proj_space_name)))

        push!(stmts, Graphviz.Edge(gp_proj1(p), Dict(:label=>"proj₁", :style=>"dashed")))
        push!(stmts, Graphviz.Edge(gp_proj2(p), Dict(:label=>"proj₂", :style=>"dashed")))

        res_name = String(d[p, :op2])
        push!(stmts, Graphviz.Edge(gp_projRes(p), Dict(:label=>res_name)))
    end
    #Need to add user access for more customizability later
    Graphviz.Graph("G", true, "neato", stmts, Dict(), Dict(:shape=>"oval"), Dict())
end


Graphics.to_graphviz(F::AbstractDecapode; kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(F; kw...))


function Catlab.Graphics.to_graphviz_property_graph(d::AbstractNamedDecapode; kw...)
    pg = PropertyGraph{Any}(;kw...)
    vids = map(parts(d, :Var)) do v
      add_vertex!(pg, label=varname(d,v))
    end

    map(parts(d, :Op1)) do op
      s, t = d[op, :src], d[op, :tgt]
      add_edge!(pg, vids[s],vids[t], label=String(d[op,:op1]))
    end

    map(parts(d, :Op2)) do op
      s, t = d[op, :proj1], d[op, :proj2]
      r = d[op, :res]
      v = add_vertex!(pg, label="$(spacename(d, s))×$(spacename(d,t))", shape="rectangle")
      add_edge!(pg, v, vids[s], label="π₁", style="dashed")
      add_edge!(pg, v, vids[t], label="π₂", style="dashed")
      add_edge!(pg, v, vids[r], label=String(d[op, :op2]))
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
        matches = findvid(G, d, u)
        length(matches) == 1 || error("did not find a unique vertex match for Σ$s")
        uG = first(first(matches))
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