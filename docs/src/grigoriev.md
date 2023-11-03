# Halfar's model of glacial flow

Let's model glacial flow using a model of how ice height of a glacial sheet changes over time, from P. Halfar's 1981 paper: "On the dynamics of the ice sheets".

Let's run the Halfar shallow ice/ shallow slope model on some "real world" data for ice thickness. Van Tricht et al. in their 2023 communication [Measuring and modelling the ice thickness of the Grigoriev ice cap (Kyrgyzstan) and comparison with global dataset](https://tc.copernicus.org/articles/17/4315/2023/tc-17-4315-2023.html) published ice thickness data on an ice cap and stored their data in a TIF. In this document, we will demonstrate how to parse such data and execute a Decapodes model on these initial conditions.

For the parameters to Glen's law, we will use those used in the [Community Ice Sheet Model benchmark](https://cise.ufl.edu/~luke.morris/cism.html). Of course, the parameters of this Kyrgyzstani ice cap likely differ from these by quite some amount, but they are a good place to start. Further, this ice cap does not satisfy the "shallow slope" assumption across the entire domain.

``` @example DEC
# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes

# External Dependencies
using FileIO  
using Interpolations
using MLStyle
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using SparseArrays
using CairoMakie
using GeometryBasics: Point2
Point2D = Point2{Float64}
Point3D = Point3{Float64}; # hide
```

# Loading a Scientific Dataset
The ice thickness data is [stored in a TIF](https://zenodo.org/api/records/7735970/files-archive). We have downloaded it locally, and load it using basic `FileIO`.

``` @example DEC
file_name = "Icethickness_Grigoriev_ice_cap_2021.tif"
ice_thickness_tif = load(file_name)
```

This data may appear to be a simple binary mask, but that is only because values with no ice are set to `-Inf`. We will account for this we interpolate our data.

We use the `Interpolations.jl` library to interpolate this dataset:

``` @example DEC
# Taking the coordinates to be from the extrema of the measured points:
const MIN_X = 243504.5
const MAX_X = 245599.8
const MIN_Y = 4648894.5
const MAX_Y = 4652179.7
ice_coords = (range(MIN_X, MAX_X, length=size(ice_thickness_tif,1)),
              range(MIN_Y, MAX_Y, length=size(ice_thickness_tif,2)))
# Note that the tif is set to -floatmax(Float32) where there is no ice.
# For our purposes, this is equivalent to 0.0.
ice_interp = LinearInterpolation(ice_coords, Float32.(ice_thickness_tif))
```

To use this interpolating object `ice_interp`, we can simply query it for the value at some coordinates: `ice_interp(x,y)`.

Let's generate a triangulated grid located at the appropriate coordinates:

``` @example DEC
include("../../examples/grid_meshes.jl")
# Specify a resolution:
RES_X = (MAX_X-MIN_X)/30.0
RES_Y = RES_X
# Generate the mesh with appropriate dimensions and resolution:
s′ = triangulated_grid(
                       MAX_X-MIN_X, MAX_Y-MIN_Y,
                       RES_X, RES_Y, Point3D)
# Shift it into place:
s′[:point] = map(x -> x + Point3D(MIN_X, MIN_Y, 0), s′[:point])
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())
wireframe(s)
```

The coordinates of a vertex are stored in `s[:point]`. Let's use our interpolator to assign ice thickness values to each vertex in the mesh:

``` @example DEC
# These are the values used by the CISM benchmark:
n = 3
ρ = 910
g = 9.8101
A = fill(1e-16, ne(s))

h₀ = map(s[:point]) do (x,y,_)
  tif_val = ice_interp(x,y)
  # Accommodate for the -∞'s that encode "no ice".
  tif_val < 0.0 ? 0.0 : tif_val
end

# Store these values to be passed to the solver.
u₀ = ComponentArray(h=h₀, stress_A=A)
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g)
```

# Defining and Composing Models
For exposition on this Halfar Decapode, see our [Glacial Flow](https://algebraicjulia.github.io/Decapodes.jl/dev/ice_dynamics) docs page. You can skip ahead to the next section.

``` @example DEC
halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯(d(h)))^(n-1)) * avg₀₁(h^(n+2)))
end

glens_law = @decapode begin
  Γ::Form1
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

ice_dynamics_composition_diagram = @relation () begin
  dynamics(h,Γ,n)
  stress(Γ,n)
end

ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:h,:Γ,:n]),
  Open(glens_law, [:Γ,:n])])

ice_dynamics = apex(ice_dynamics_cospan)
to_graphviz(ice_dynamics)
```

# Define our functions

``` @example DEC
include("sharp_op.jl")
function generate(sd, my_symbol; hodge=GeometricHodge())
  ♯_m = ♯_mat(sd, AltPPSharp())
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      append!(J, [s[e,:∂v0],s[e,:∂v1]])
      append!(I, [e,e])
      append!(V, [0.5, 0.5])
  end
  avg_mat = sparse(I,J,V)
  op = @match my_symbol begin
    :♯ => x -> begin
      ♯(sd, EForm(x))
    end
    :mag => x -> begin
      norm.(x)
    end
    :avg₀₁ => x -> begin
      avg_mat * x
    end
    :^ => (x,y) -> x .^ y
    :* => (x,y) -> x .* y
    :abs => x -> abs.(x)
    :show => x -> begin
      println(x)
      x
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

# Generate simulation

``` @example DEC
sim = eval(gensim(ice_dynamics2D, dimension=2))
fₘ = sim(s, generate)
```

# Run

``` @example DEC
tₑ = 1e1

@info("Solving Grigoriev Ice Cap")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
@save "grigoriev.jld2" soln
```

# Visualize

``` @example DEC
# Visualize the initial conditions.
function plot_ic()
  f = Figure()
  ax = Axis(f[1,1],
            title="Grigoriev Ice Cap Initial Thickness [m]",
            xticks = range(MIN_X, MAX_X; length=5),
            yticks = range(MIN_Y, MAX_Y; length=5))
  msh = mesh!(ax, s′, color=soln(0.0).h, colormap=:jet)
  Colorbar(f[1,2], msh)
  f
end
f = plot_ic()
save("grigoriev_ic.png", f)

# Visualize the final conditions.
function plot_fc()
  f = Figure()
  ax = Axis(f[1,1],
            title="Grigoriev Ice Cap Final Thickness [m]",
            xticks = range(MIN_X, MAX_X; length=5),
            yticks = range(MIN_Y, MAX_Y; length=5))
  msh = mesh!(ax, s′, color=soln(tₑ).h, colormap=:jet)
  Colorbar(f[1,2], msh)
  f
end
f = plot_fc()
save("grigoriev_fc.png", f)

# Create a gif
function save_dynamics(save_file_name)
  time = Observable(0.0)
  h = @lift(soln($time).h)
  f,a,o = mesh(s′, color=h, colormap=:jet,
             colorrange=extrema(soln(tₑ).h);
             axis = (; title = @lift("Grigoriev Ice Cap Dynamic Thickness [m] at time $($time))")))
  Colorbar(f[1,2], limits=extrema(soln(0.0).h), colormap=:jet)
  timestamps = range(0, tₑ, step=1e-1)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics("grigoriev.gif")
```

![Grigoriev_ICs](grigoriev_ic.png)
![Grigoriev_FCs](grigoriev_fc.png)
![Grigoriev_Dynamics](grigoriev.gif)
