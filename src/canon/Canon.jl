module Canon

using DiagrammaticEquations
using DiagrammaticEquations.Deca

using Markdown

export @docapode

function create_pode_expr(t) 
  modelname, source, desc, variable, pode = t
  variable = Symbol(variable)
  Base.remove_linenums!(pode)

  modeldef = sprint.(Base.show_unquoted,pode.args)
  modeldef = join(map(modeldef) do line; "    $line
                  " end, "\n")

  docstring = Markdown.parse("""

  # $modelname

  [Source]($source)

  $desc

  ## Model 

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

end
