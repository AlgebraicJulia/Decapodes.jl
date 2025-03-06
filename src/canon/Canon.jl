module Canon

using DiagrammaticEquations
using DiagrammaticEquations.Deca

using DocumenterCitations
using Markdown

bib = CitationBibliography(
    joinpath(@__DIR__, "../../docs", "src", "decapodes_documenter.bib");
    style=:numeric)


struct DecapodeModel
    value::Symbol
end

export @docapode

function create_pode_expr(t) 
  modelname, source, desc, variable, pode = t
  variable = Symbol(variable)
 
  modeldef = String[]; lnns = LineNumberNode[]
  foreach(pode.args) do arg
      if arg isa LineNumberNode
          !isempty(lnns) && arg.line - last(lnns).line > 1 ? push!(modeldef, "") : nothing
          push!(lnns, arg)
      else
          push!(modeldef, sprint(Base.show_unquoted, arg))
      end
  end
  modeldef = join(map(modeldef) do line; "    $line" end, "\n")

  docstring = Markdown.parse("""

  **$modelname**

  $(!isempty(source) ? """[Source]($source)
    """ : "")
  $desc

  **Model** 

  $modeldef

  """)
  quote 
    export $variable
    $variable = @decapode $pode
    @doc $docstring $variable
  end
end


macro docapode(name, source, description, variable, pode)  
    create_pode_expr((name, source, description, variable, pode)) |> esc
end

include("Physics.jl")
include("Chemistry.jl")
include("Biology.jl")
include("Environment.jl")
include("Oncology.jl")

end
