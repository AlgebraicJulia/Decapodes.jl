# Couple Ice and Water Dynamics

Let's use Decapodes to implement the incompressible Navier-Stokes as given by [Mohamed et al.](https://arxiv.org/abs/1508.01166). We will run these dynamics [on the sphere](https://youtu.be/k0hFhAvhHvs?si=Wi9-OgBbAODtxMtb). We will couple this model with Halfar glacier dynamics [on the sphere](https://algebraicjulia.github.io/Decapodes.jl/dev/ice_dynamics/#2-Manifold-in-3D). For the initial conditions of the Halfar ice thickness, we will use scientific dataset for Greenland, much like the scientific dataset used for the [Grigoriev ice cap](https://algebraicjulia.github.io/Decapodes.jl/dev/grigoriev/) simulation.

Note that the time scale at which ice creeps is much larger than the time scale at which the water in the ocean would flow. So we can either choose to model a very slow moving fluid around the ice (like a storm on a gas giant), or we can choose to model on a shorter timescale, on which the ice does not move very much.

```@example DEC
# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics: Point3
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
Point3D = Point3{Float64};
```

## Specify our models.

Our first component is the Mohamed et al. formulation of the incompressible Navier-Stokes equations. We will call the flow here "w". This will be the flow after collisions with glaciers are considered.

This is [Equation 10](https://arxiv.org/abs/1508.01166) for N=2.
```@example DEC
eq10forN2 = @decapode begin
  (𝐮,w)::DualForm1
  (P, 𝑝ᵈ)::DualForm0
  μ::Constant

  𝑝ᵈ == P + 0.5 * ι₁₁(w,w)

  ∂ₜ(𝐮) == μ * ∘(d, ⋆, d, ⋆)(w) + (-1)*⋆₁⁻¹(∧ᵈᵖ₁₀(w, ⋆(d(w)))) + d(𝑝ᵈ)
end
to_graphviz(eq10forN2)
```

Halfar's equation and Glen's law are composed like so:

```@example DEC
halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ∂ₜ(h) == ∘(⋆, d, ⋆)(Γ  * d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end

glens_law = @decapode begin
  Γ::Form1
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

ice_dynamics_composition_diagram = @relation () begin
  dynamics(Γ,n)
  stress(Γ,n)
end

ice_dynamics = apex(oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:Γ,:n]),
   Open(glens_law, [:Γ,:n])]))

to_graphviz(ice_dynamics, verbose=false)
```

We now have our dynamics that govern glaciers, and our dynamics that govern water. We need to specify the physics of what happens when glaciers and water interact.
There are many options, and the choice you make depends on the time-scale and resolution of the dynamics that you are interested in.

An interaction between glacier and water dynamics can look like this:

```@example DEC
ice_water_composition_diagram = @relation () begin
  glacier_dynamics(ice_thickness)
  water_dynamics(flow, flow_after)

  interaction(ice_thickness, flow, flow_after)
end
to_graphviz(ice_water_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")
```
, where `flow_after` is the flow of water after interaction with ice is considered.

We will use the language of Decapodes to encode the dynamics that ice blocks water from flowing.

We can detect the ice with a sigmoid function. Where there is ice, we want the flow to be 0, and where there is no ice, we will not impede the flow. We won't consider any further special boundary conditions between ice and water here. Since h is a scalar-like quantity, and flow is a vector-like quantity, we can relate them using the wedge product operator from the exterior calculus. We can state these dynamics using the language of Decapodes like so:

```@example DEC
blocking = @decapode begin
  h::Form0
  (𝐮,w)::DualForm1

  w == (1-σ(h)) ∧ᵖᵈ₀₁ 𝐮
end
to_graphviz(blocking)
```
, where σ is a sigmoid function that is 0 when d(h) is 0, and goes to 1 otherwise. We see that w is indeed defined as 𝐮, after interacting with the ice boundary is considered.

We can apply our composition diagram to generate our physics:
```@example DEC
ice_water = apex(oapply(ice_water_composition_diagram,
  [Open(ice_dynamics, [:dynamics_h]),
   Open(eq10forN2,    [:𝐮, :w]),
   Open(blocking,     [:h, :𝐮, :w])]))

to_graphviz(ice_dynamics, verbose=false)
```

We can now generate our simulation:

```@example DEC
open("ice_water.jl", "w") do f
  write(f, string(gensim(expand_operators(ice_water))))
end
sim = include("ice_water.jl")
```

## Meshes and Initial Conditions

Since we want to demonstrate these physics on the Earth, we will use one of our icosphere discretizations with the appropriate radius.

``` @example DEC
rₑ = 6378e3 # [km]
s′ = loadmesh(Icosphere(5, rₑ))
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())
wireframe(s)
```

Let's demonstrate how to add operators by providing the definition of a sigmoid function:

```@example DEC
sigmoid(x) = (2 ./ (1 .+ exp.(-x*1e2)) .- 1)
function generate(s, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    # This is a new function.
    :σ => sigmoid
    :mag => x -> norm.(x)
    :^ => (x,y) -> x .^ y
    # Remaining operations (such as our differential operators) are built-in.
    _ => default_dec_matrix_generate(s, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end;
```

Let's combine our mesh with our physics to instantiate our simulation:

```@example DEC
fₘ = sim(s, generate);
```

We can now supply initial conditions:

```@example DEC
#TODO: Grab ice data from interpolator.

ice_thickness = map(s[:point]) do (_,_,z)
  z < 0.8*rₑ ? 0 : 1
end

flow = dec_dual_derivative(0,s) *
  map(s[s[:tri_center], :dual_point]) do (_,_,z)
    (rₑ-abs(z))/rₑ
  end
# There is no water "under" the ice:
flow = dec_wedge_product_pd(Tuple{0,1},s)(1 .- sigmoid(ice_thickness), flow)

u₀ = ComponentArray(
  ice_thickness = ice_thickness,
  flow = flow,
  water_dynamics_P = zeros(ntriangles(s)))

constants_and_parameters = (
  glacier_dynamics_n = 3,
  glacier_dynamics_stress_A = fill(1e-16, ne(s)),
  glacier_dynamics_stress_ρ = 910,
  glacier_dynamics_stress_g = 9.8101,
  water_dynamics_μ = 0.01);
```
  
## Execute the Simulation

We specified our physics, our mesh, and our initial conditions. We have everything we need to execute the simulation.

```@example DEC
tₑ = 1e2

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Vern7())
soln.retcode != :Unstable || error("Solver was not stable")

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Vern7())
@show soln.retcode
@info("Done")
```

## Results

In the DEC, vorticity is encoded with `d⋆`, and speed can be encoded with `norm ♯`. We can use our operators from CombinatorialSpaces.jl to create our GIFs.

```@example DEC
ihs0 = dec_inv_hodge_star(Val{0}, s, GeometricHodge())
dd1 = dec_dual_derivative(1, s)
♯_m = ♯_mat(s, LLSDDSharp())
using LinearAlgebra: norm
function vorticity(α)
  ihs0*dd1*α
end
function speed(α)
  norm.(♯_m * α)
end

function save_vorticity(is_2d=false)
  frames = 200
  time = Observable(0.0)
  fig = Figure(title = @lift("Vorticity at $($time)"))
  ax = is_2d ?
    CairoMakie.Axis(fig[1,1]) :
    LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.mesh!(ax, s′,
    color=@lift(vorticity(soln($time).flow)),
    colorrange=extrema(vorticity(soln(tₑ).flow)).*.9,
    colormap=:jet)

  Colorbar(fig[1,2], msh)
  record(fig, "vorticity_ice_water.gif", range(0.0, tₑ; length=frames); framerate = 20) do t
    time[] = t
  end
end
save_vorticity(false)
```

![Vorticity](vorticity_ice_water.gif)

Let's look at the dynamics of the ice as well:

``` @example DEC
begin
  frames = 200
  fig = Figure()
  ax = LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.mesh!(ax, s′, color=soln(0).ice_thickness, colormap=:jet, colorrange=extrema(soln(0).ice_thickness))

  Colorbar(fig[1,2], msh)
  record(fig, "halmo_ice.gif", range(0.0, tₑ; length=frames); framerate = 20) do t
    msh.color = soln(t).ice_thickness
  end
end
```

![HalfarMohamedIce](halmo_ice.gif)

