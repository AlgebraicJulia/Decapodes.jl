module OpenDiagrams
export OpenDiagram, oapply

using Catlab, Catlab.CategoricalAlgebra, Catlab.WiringDiagrams, Catlab.Programs
import Catlab.CategoricalAlgebra: apex, legs, feet
import Catlab.WiringDiagrams: oapply
using Catlab.Programs.DiagrammaticPrograms: NamedGraph # FIXME: Should export?

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

end
