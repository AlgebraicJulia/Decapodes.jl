```@meta
EditURL = "../../../literate/oncology/tumor_proliferation_invasion.jl"
```

Load Dependencies

````@example tumor_proliferation_invasion
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations, DiagrammaticEquations.Deca
using Distributions
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using CairoMakie
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}
````

Load in our Decapodes models

````@example tumor_proliferation_invasion
using Decapodes.Canon.Oncology
````

Let's examine our models. Here's the tumor invasion model:

````@example tumor_proliferation_invasion
to_graphviz(invasion)
````

The Logistic growth model

````@example tumor_proliferation_invasion
to_graphviz(logistic)
````

The Gompertz growth model

````@example tumor_proliferation_invasion
to_graphviz(gompertz)
````

Load in a mesh, initial conditions, and a plotting function

````@example tumor_proliferation_invasion
function show_heatmap(Cdata)
  heatmap(reshape(Cdata, (floor(Int64, sqrt(length(Cdata))), floor(Int64, sqrt(length(Cdata))))))
end

s = triangulated_grid(50,50,0.2,0.2,Point2D);
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point2D}(s);
subdivide_duals!(sd, Circumcenter());

constants_and_parameters = (
  invasion_Dif = 0.005,
  invasion_Kd = 0.5,
  Cmax = 10)
````

"The model ... considers an equivalent radially symmetric tumour"
- Murray J.D., Glioblastoma brain tumours

````@example tumor_proliferation_invasion
c_dist  = MvNormal([25, 25], 2)
C = 100 * [pdf(c_dist, [p[1], p[2]]) for p in sd[:point]]

u₀ = ComponentArray(C=C)
````

Compose our Proliferation-Invasion models

````@example tumor_proliferation_invasion
proliferation_invasion_composition_diagram = @relation () begin
  proliferation(C, fC, Cmax)
  invasion(C, fC, Cmax)
end

logistic_proliferation_invasion_cospan = oapply(proliferation_invasion_composition_diagram,
  [Open(logistic, [:C, :fC, :Cmax]),
   Open(invasion, [:C, :fC, :Cmax])])

logistic_proliferation_invasion = apex(logistic_proliferation_invasion_cospan)

gompertz_proliferation_invasion_cospan = oapply(proliferation_invasion_composition_diagram,
  [Open(gompertz, [:C, :fC, :Cmax]),
   Open(invasion, [:C, :fC, :Cmax])])

gompertz_proliferation_invasion = apex(gompertz_proliferation_invasion_cospan)
````

Generate the logistic simulation

````@example tumor_proliferation_invasion
logistic_sim = evalsim(logistic_proliferation_invasion)

lₘ = logistic_sim(sd, default_dec_generate, DiagonalHodge())
````

Execute the logistic simulation

````@example tumor_proliferation_invasion
tₑ = 15.0

prob = ODEProblem(lₘ, u₀, (0, tₑ), constants_and_parameters)
logistic_soln = solve(prob, Tsit5())

show_heatmap(logistic_soln(tₑ).C)
````

Generate the Gompertz simulation

````@example tumor_proliferation_invasion
gompertz_sim = evalsim(gompertz_proliferation_invasion)
gₘ = gompertz_sim(sd, default_dec_generate, DiagonalHodge())
````

Execute the Gompertz simulation

````@example tumor_proliferation_invasion
prob = ODEProblem(gₘ, u₀, (0, tₑ), constants_and_parameters)
gompertz_soln = solve(prob, Tsit5())

show_heatmap(gompertz_soln(tₑ).C)
````

