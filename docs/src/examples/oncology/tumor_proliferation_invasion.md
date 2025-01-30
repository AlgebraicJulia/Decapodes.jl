```@meta
EditURL = "../../../literate/oncology/tumor_proliferation_invasion.jl"
```

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

Let's examine our models. Here's the tumor invasion model with the logistic and Gompertz growth models.

````@example tumor_proliferation_invasion
@doc invasion
````

````@example tumor_proliferation_invasion
@doc logistic
````

````@example tumor_proliferation_invasion
@doc gompertz
````

Load in a mesh and a plotting function

````@example tumor_proliferation_invasion
function show_heatmap(Cdata)
  heatmap(reshape(Cdata, (floor(Int64, sqrt(length(Cdata))), floor(Int64, sqrt(length(Cdata))))))
end

mesh = triangulated_grid(50,50,0.2,0.2,Point2D);
dualmesh = EmbeddedDeltaDualComplex2D{Bool, Float64, Point2D}(mesh);
subdivide_duals!(dualmesh, Circumcenter());
nothing #hide
````

Let's define initial conditions

````@example tumor_proliferation_invasion
constants_and_parameters = (invasion_Dif = 0.005, invasion_Kd = 0.5, Cmax = 10)
````

Here we follow the assumption "The model ... considers an equivalent radially symmetric tumour", Murray J.D., Glioblastoma brain tumours, by initializing the tumor with a normal distribution.

````@example tumor_proliferation_invasion
c_dist  = MvNormal([25, 25], 2)
C = 100 * [pdf(c_dist, [p[1], p[2]]) for p in dualmesh[:point]]
u₀ = ComponentArray(C=C)
````

Let's define how our Proliferation-Invasion models will relate to one another.

````@example tumor_proliferation_invasion
proliferation_invasion_composition_diagram = @relation () begin
  proliferation(C, fC, Cmax)
  invasion(C, fC, Cmax)
end
````

Now let's specify which sub-models slot into our system. We use the same pattern for two different models: the first model pertains to a logistic growth model,

````@example tumor_proliferation_invasion
logistic_proliferation_invasion_cospan = oapply(proliferation_invasion_composition_diagram,
  [Open(logistic, [:C, :fC, :Cmax]),
   Open(invasion, [:C, :fC, :Cmax])])
logistic_proliferation_invasion = apex(logistic_proliferation_invasion_cospan)
````

The second model uses the same composition pattern but swaps out the logistic growth mode for a Gompertz growth model.

````@example tumor_proliferation_invasion
gompertz_proliferation_invasion_cospan = oapply(proliferation_invasion_composition_diagram,
  [Open(gompertz, [:C, :fC, :Cmax]),
   Open(invasion, [:C, :fC, :Cmax])])
gompertz_proliferation_invasion = apex(gompertz_proliferation_invasion_cospan)
````

Generate the logistic simulation

````@example tumor_proliferation_invasion
logistic_sim = evalsim(logistic_proliferation_invasion)
lₘ = logistic_sim(dualmesh, default_dec_generate, DiagonalHodge())
````

Execute the logistic simulation

````@example tumor_proliferation_invasion
tₑ = 15.0
problem = ODEProblem(lₘ, u₀, (0, tₑ), constants_and_parameters)
logistic_solution = solve(problem, Tsit5());
nothing #hide
````

Let's examine the solution using the heatmap equation we defined.

````@example tumor_proliferation_invasion
show_heatmap(logistic_solution(tₑ).C)
````

Generate the Gompertz simulation

````@example tumor_proliferation_invasion
gompertz_sim = evalsim(gompertz_proliferation_invasion)
gₘ = gompertz_sim(dualmesh, default_dec_generate, DiagonalHodge())
````

Execute the Gompertz simulation

````@example tumor_proliferation_invasion
problem = ODEProblem(gₘ, u₀, (0, tₑ), constants_and_parameters)
gompertz_solution = solve(problem, Tsit5());
nothing #hide
````

Let's examine this solution now.

````@example tumor_proliferation_invasion
show_heatmap(gompertz_solution(tₑ).C)
````

