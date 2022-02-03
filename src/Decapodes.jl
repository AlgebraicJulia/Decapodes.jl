module Decapodes

using Requires

include("Diagrams.jl")
include("OpenDiagrams.jl")
include("Schedules.jl")
include("Simulations.jl")
include("Examples.jl")

function __init__()
  @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("Debug.jl")
end
end
