module Canon

using DiagrammaticEquations
using DiagrammaticEquations.Deca

using DocumenterCitations
using Markdown

bib = CitationBibliography(
    joinpath(@__DIR__, "../../docs", "src", "decapodes_documenter.bib");
    style=:numeric)

"""
    format_citation_key(entries, key)

Format a single citation key as "Author (Year)" using the bibliography entries.
Falls back to the raw key if it is not found.
"""
function format_citation_key(entries, key)
  haskey(entries, key) || return key
  entry = entries[key]
  author = if !isempty(entry.authors)
    entry.authors[1].last
  else
    key
  end
  year = entry.date.year
  isempty(year) ? author : "$author ($year)"
end

"""
    expand_citations(desc, entries)

Replace DocumenterCitations `[key1, key2](@cite)` patterns in `desc` with
formatted "Author (Year)" text so that docstrings render properly at runtime.
"""
function expand_citations(desc, entries)
  cite_re = r"\[([^\]]+)\]\(@cite\)"
  replace(desc, cite_re => function(m)
    # Strip surrounding "[" and "](@cite)" to get the inner key(s)
    inner = m[2:end-8]  # remove leading '[' and trailing '](@cite)'
    keys = strip.(split(inner, ","))
    join([format_citation_key(entries, k) for k in keys], ", ")
  end)
end

struct DecapodeModel
    value::Symbol
end

export @docapode

function create_pode_expr(t) 
  modelname, source, desc, variable, pode = t
  desc = expand_citations(desc, bib.entries)
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
