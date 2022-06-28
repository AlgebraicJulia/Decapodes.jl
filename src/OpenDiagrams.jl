module OpenDiagrams
export OpenDiagram, oapply, draw_diagram

using Catlab, Catlab.CategoricalAlgebra, Catlab.WiringDiagrams, Catlab.Programs
import Catlab.CategoricalAlgebra: apex, legs, feet
import Catlab.WiringDiagrams: oapply
using Catlab.Programs.DiagrammaticPrograms: NamedGraph # FIXME: Should export?
using Catlab.Graphics

using CombinatorialSpaces.ExteriorCalculus: ExtCalc2D
using CombinatorialSpaces
using Catlab.CategoricalAlgebra.Categories: Functor, Cat
using Catlab.CategoricalAlgebra.FinCats: FinCatSize, FinCatPresentation

#const Diagram2D = FinFunctor{Dom, Codom} where {Ob, Hom, Dom, Codom <: FinCatPresentation{ExtCalc2D, Ob, Hom}}

const Diagram2D = Functor{Dom, Codom} where
                    {Ob, Hom, Dom<:(Cat{Ob, Hom, FinCatSize} where
                      {Ob, Hom}), Codom<:FinCatPresentation{ExtCalc2D, Ob, Hom}}

""" Open diagram as a structured multicospan in R-form.

An open diagram is a diagram, represented as a `FinDomFunctor`, together with a
cospan of finite sets whose apex is the object set of the diagram's indexing
category.
"""
struct OpenDiagram{F<:FinDomFunctor,Cosp<:Multicospan}
  functor::F
  cospan::Cosp
end

apex(diagram::OpenDiagram) = diagram.functor
legs(diagram::OpenDiagram) = legs(diagram.cospan)
feet(diagram::OpenDiagram) = feet(diagram.cospan)

function OpenDiagram(F::FinDomFunctor, names::AbstractVector{Symbol})
  g = graph(dom(F))
  legs = map(name -> FinFunction([incident(g, name, :vname)], nv(g)), names)
  OpenDiagram(F, Multicospan(FinSet(1), legs))
end

const OpenFreeDiagramOb, OpenFreeDiagram = OpenACSetTypes(FreeDiagram, :V)

OpenFreeDiagram(d::FreeDiagram{Ob,Hom}, args...) where {Ob,Hom} =
  OpenFreeDiagram{Ob,Hom}(d, args...)

function oapply(uwd::RelationDiagram, diagrams::AbstractVector{<:OpenDiagram})
  open_free_diagram = oapply(uwd, map(diagrams) do d
    OpenFreeDiagram(FreeDiagram(apex(d)), d.cospan)
  end)
  free_diagram = apex(open_free_diagram)
  leg_functions = map(leg -> leg[:V], legs(open_free_diagram))

  g = NamedGraph{Symbol,Union{Symbol,Nothing}}()
  copy_parts!(g, free_diagram)
  g[:vname] = collect(apex(colimit_vertex_names(uwd, diagrams)))
  g[:ename] = nothing # FIXME: Coproduct edge names, allowing for nulls.

  diagram = FinFunctor(ob(free_diagram), hom(free_diagram),
                       FinCat(g), codom(apex(first(diagrams))))
  OpenDiagram(diagram, Multicospan(leg_functions))
end

function colimit_vertex_names(uwd::RelationDiagram,
                              diagrams::AbstractVector{<:OpenDiagram})
  @assert nboxes(uwd) == length(diagrams)
  bfd = BipartiteFreeDiagram{FinSet{<:Any,Symbol},FinFunction}()
  add_vertices₁!(bfd, njunctions(uwd), ob₁=map(uwd[:variable]) do var
    FinSet([var])
  end)
  add_vertices₂!(bfd, nboxes(uwd), ob₂=map(diagrams) do diagram
    FinSet(graph(dom(diagram.functor))[:vname])
  end)
  for (b, diagram) in zip(boxes(uwd), diagrams)
    g = graph(dom(diagram.functor))
    for (p, leg) in zip(ports(uwd, b), legs(diagram.cospan))
      j, v = junction(uwd, p), only(collect(leg))
      f = FinFunction(Dict(uwd[j, :variable] => g[v, :vname]), ob₂(bfd, b))
      add_edge!(bfd, j, b, hom=f)
    end
  end
  colimit(bfd)
end

####################
# Graphics Tooling #
####################

draw_diagram(d::Diagram2D) = to_graphviz(d, node_labels=true,
                                         prog="neato",
                                         node_attrs=Dict(:shape=>"oval"),
                                         graph_attrs=Dict(:nodesep=>"4.0"))
using Catlab.Graphs
using Catlab.Graphs.BasicGraphs
function to_graph(J::FinCat)
    g = BasicGraphs.Graph()
    copy_parts!(g, graph(J))
    return g
end

Graphics.to_graphviz(F::Diagram2D; kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(F; kw...))

function GraphvizGraphs.to_graphviz_property_graph(F::Diagram2D; kw...)
    simplify_vname(s) = begin
      table = Dict(
    "Form0(X)" => "Ω₀",
    "Form1(X)" => "Ω₁",
    "Form2(X)" => "Ω₂",
    "DualForm0(X)" => "Ω̃₀",
    "DualForm1(X)" => "Ω̃₁",
    "DualForm2(X)" => "Ω̃₂",
    "otimes(Form0(X),Form0(X))" => "Ω₀²",
    "otimes(Form1(X),Form1(X))" => "Ω₁²",
    "otimes(Form1(X),DualForm2(X))" => "Ω₁×Ω̃₂",
    "otimes(Form1(X),Form1(X),Form1(X))" => "Ω₁³",
    "otimes(Form1(X),Form1(X),Form1(X),Form1(X))" => "Ω₁⁴"
    )
    if string(s) in keys(table)
        return table[string(s)]
    else
        println(string(s))
        return string(s)
    end
    end

    simplify_ename(s) = begin
        b = IOBuffer()
        show_unicode(b, s)
        return replace(String(take!(b)), r"{X}"=>"")
    end

    J = dom(F)
    G = graph(J)
    pg = GraphvizGraphs.to_graphviz_property_graph(to_graph(J); kw...)
    for v in vertices(G)
        lᵥ = G[v, :vname]
        tᵥ = simplify_vname(F.ob_map[v])
        set_vprops!(pg, v, Dict(:label => "$(lᵥ):$tᵥ"))
    end
    for e in edges(G)
        tₑ = F.hom_map[e]
        set_eprops!(pg, e, Dict(:label => "$(simplify_ename(tₑ))"))
    end
    pg
end


Graphics.to_graphviz(cosp::OpenDiagram; kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(cosp; kw...))


function GraphvizGraphs.to_graphviz_property_graph(cosp::OpenDiagram; kw...)
    pg = GraphvizGraphs.to_graphviz_property_graph(cosp.functor; kw...)
    label(I, l, i) = length(dom(l)) > 1 ? "$(i):$I" : "$I"
    for (I,l) in enumerate(legs(cosp.cospan))
        for i in dom(l)
            v = add_vertex!(pg)
            set_vprops!(pg, v, Dict(:label=>label(I, l,i), :shape=>"circle", :color=>"blue"))
            e = add_edge!(pg, v, l(i))
            set_eprops!(pg, e, Dict(:style=>"dashed"))
        end
    end
    return pg
end

end
