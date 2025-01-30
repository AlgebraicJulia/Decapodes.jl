# Halfar-EBM-Water

This docs page demonsrates a composition of the Halfar model of ice dynamics with the Budyko-Sellers energy-balance model (EBM) defining temperature dynamics. Surface temperature affects the rate at which ice diffuses, and melts ice into a water density term. This water density then diffuses across the domain.

We execute these dynamics on real sea ice thickness data provided by the Polar Science Center.

``` @example DEC
# Import Dependencies 

# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes

# External Dependencies
using ComponentArrays
using CoordRefSystems
using CairoMakie
using LinearAlgebra
using MLStyle
using NearestNeighbors
using NetCDF
using OrdinaryDiffEq
Point3D = Point3{Float64}
```

First, we load a mesh. We will execute these dynamics on the sphere:

``` @example DEC
# Load a mesh
s_plots = loadmesh(Icosphere(7));
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s_plots);
subdivide_duals!(s, Barycenter());
wireframe(s_plots)
```

## Load data

The data provided by the Polar Science Center is given as a NetCDF file. Ice thickness is a matrix with the same dimensions as a matrix provided Latitude and Longitude of the associated point on the Earth's surface. We need to convert between polar and Cartesian coordinates to use this data on our mesh.

``` @example DEC
# This data can be downloaded from source here:
# https://pscfiles.apl.uw.edu/axel/piomas20c/v1.0/monthly/piomas20c.heff.1901.2010.v1.0.nc
ice_thickness_file = "piomas20c.heff.1901.2010.v1.0.nc"
run(`curl -o $ice_thickness_file https://cise.ufl.edu/"~"luke.morris/piomas20c.heff.1901.2010.v1.0.nc`)

# Use ncinfo(ice_thickness_file) to interactively get information on variables.
# Sea ice thickness ("sit") has dimensions of [y, x, time].
# y,x index into "Latitude" and "Longitude" variables.
# Time is in units of days since 1901-01-01.
lat = ncread(ice_thickness_file, "Latitude")
lon = ncread(ice_thickness_file, "Longitude")
sit = ncread(ice_thickness_file, "sit")

# Convert latitude from [90, -90] to [0, 180] for convenience.
lat .= -lat .+ 90

# Convert mesh points from Cartesian to spherical coordinates.
p_sph = map(point(s)) do p
  p = convert(Spherical, Cartesian(p...))
  [rad2deg(p.θ).val, rad2deg(p.ϕ).val]
end

# Note: You can instead use an algebraic parameterization, rather than nearest-neighbor interpolation.
lat, lon = lat[:], lon[:]
ll = hcat(lat, lon)'
kdt = KDTree(ll)
sit_sph_idxs = map(p_sph) do p
  nn(kdt, p)[1]
end

sit_sph = map(sit_sph_idxs, p_sph) do i, p
  ((p[1] > maximum(lat)) || isnan(sit[i])) ? 0.0f0 : sit[i]
end

f = Figure()
ax = LScene(f[1,1], scenekw=(lights=[],))
update_cam!(ax.scene, Vec3f(0,0,0.8), Vec3f(0,0,0), Vec3f(0, 1, 1))
msh = mesh!(ax, s_plots, color=sit_sph, colormap=Reverse(:redsblues))
Colorbar(f[1,2], msh)
f
```

## Define the model

Here, the Halfar model of ice dynamics is recalled, as well as the Budyko-Sellers EBM. These these models are composed individually. They are then coupled together via `warming` and `melting` components.

``` @example DEC
halfar_eq2 = @decapode begin
  (h, melt)::Form0
  Γ::Form1
  n::Constant

  ∂ₜ(h)  == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯(d(h)))^(n-1)) * avg₀₁(h^(n+2))) - melt
end

glens_law = @decapode begin
  (A,Γ)::Form1
  (ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

ice_dynamics_composition_diagram = @relation () begin
  dynamics(Γ,n)
  stress(Γ,n)
end

ice_dynamics = apex(oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:Γ,:n]),
   Open(glens_law, [:Γ,:n])]))

draw_composition(ice_dynamics_composition_diagram)
```

The composition pattern tells you how to couple the variables and introduces namespaces that we will use later when supplying initial conditions.

The following code creates the Budyko-Sellers model as a composite of individual terms.

```@example DEC

energy_balance = @decapode begin
  (Tₛ, ASR, OLR, HT)::Form0
  C::Constant

  ∂ₜ(Tₛ) == (ASR - OLR + HT) ./ C
end

absorbed_shortwave_radiation = @decapode begin
  (Q, ASR)::Form0
  α::Constant

  ASR == (1 .- α) .* Q
end

outgoing_longwave_radiation = @decapode begin
  (Tₛ, OLR)::Form0
  (A,B)::Constant

  OLR == A .+ (B .* Tₛ)
end

heat_transfer = @decapode begin
  (HT, Tₛ)::Form0
  (D,cosϕᵖ,cosϕᵈ)::Constant

  HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
end

insolation = @decapode begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end

budyko_sellers_composition_diagram = @relation () begin
  energy(Tₛ, ASR, OLR, HT)
  absorbed_radiation(Q, ASR)
  outgoing_radiation(Tₛ, OLR)
  diffusion(Tₛ, HT, cosϕᵖ)
  insolation(Q, cosϕᵖ)
end

budyko_sellers = apex(oapply(budyko_sellers_composition_diagram,
  [Open(energy_balance, [:Tₛ, :ASR, :OLR, :HT]),
   Open(absorbed_shortwave_radiation, [:Q, :ASR]),
   Open(outgoing_longwave_radiation, [:Tₛ, :OLR]),
   Open(heat_transfer, [:Tₛ, :HT, :cosϕᵖ]),
   Open(insolation, [:Q, :cosϕᵖ])]))

draw_composition(budyko_sellers_composition_diagram)
```

Our full model can then be composed by adding terms for melting of water. We will assume that the meltwater is transported by diffusion because the transport of meltwater is so much faster than the melting process itself. If you wanted to increase the physical realism of this model, using a different model of melting and water transport would be a good place to start.

``` @example DEC
warming = @decapode begin
  Tₛ::Form0
  A::Form1

  A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.01e-19)
end

melting = @decapode begin
  (Tₛ, h, melt, water)::Form0
  Dₕ₂ₒ::Constant

  melt == (Tₛ - 15)*1e-16*h
  ∂ₜ(water) == melt + Dₕ₂ₒ*Δ(water)
end

budyko_sellers_halfar_water_composition_diagram = @relation () begin
  budyko_sellers(Tₛ)
  warming(A, Tₛ)
  melting(Tₛ, h, melt)
  halfar(A, h, melt)
end

draw_composition(budyko_sellers_halfar_water_composition_diagram)
```

``` @example DEC
budyko_sellers_halfar_water = apex(oapply(budyko_sellers_halfar_water_composition_diagram,
  [Open(budyko_sellers, [:Tₛ]),
   Open(warming, [:A, :Tₛ]),
   Open(melting, [:Tₛ, :h, :melt]),
   Open(ice_dynamics, [:stress_A, :dynamics_h, :dynamics_melt])]))
nothing # hide
```

## Define initial conditions

The initial data must be specified for state variables, as well as constants and parameters.

``` @example DEC
# This is a primal 0-form, with values at vertices.
cosϕᵖ = map(x -> cos(x[1]), point(s))
# This is a dual 0-form, with values at edge centers.
cosϕᵈ = map(edges(s)) do e
  (cos(point(s, src(s, e))[1]) + cos(point(s, tgt(s, e))[1])) / 2
end

α₀ = 0.354
α₂ = 0.25
α = map(point(s)) do ϕ
  α₀ + α₂*((1/2)*(3*ϕ[1]^2 - 1))
end
A = 210
B = 2
f = 0.70
ρ = 1025
cw = 4186
H = 70
C = map(point(s)) do ϕ
  f * ρ * cw * H
end
D = 0.6

# Isothermal initial conditions:
Tₛ₀ = map(point(s)) do ϕ
  15.0
end

water = map(point(s)) do _
  0.0
end

Dₕ₂ₒ = 1e-16

n = 3
halfar_ρ = 910
g = 9.8

h₀ = sit_sph
# Store these values to be passed to the solver.
u₀ = ComponentArray(
  Tₛ = Tₛ₀,
  h = h₀,
  melting_water = water)

# The underscore-separated words are the namespaces that were introduced by oapply.
constants_and_parameters = (
  budyko_sellers_absorbed_radiation_α = α,
  budyko_sellers_outgoing_radiation_A = A,
  budyko_sellers_outgoing_radiation_B = B,
  budyko_sellers_energy_C = C,
  budyko_sellers_diffusion_D = D,
  budyko_sellers_cosϕᵖ = cosϕᵖ,
  budyko_sellers_diffusion_cosϕᵈ = cosϕᵈ,
  halfar_n = n,
  halfar_stress_ρ = halfar_ρ,
  halfar_stress_g = g,
  melting_Dₕ₂ₒ = Dₕ₂ₒ)

# Define how symbols map to Julia functions
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :♯ => begin
      sharp_mat = ♯_mat(sd, AltPPSharp())
      x -> sharp_mat * x
    end
    :mag => x -> begin
      norm.(x)
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

## Generate and run simulation 

The composed model is generated and executed for 100 years. The model is run twice to demonstrate the speed of the model after the simulation code is precompiled by the first run.

We will save the final ice thickness data in a .jld2 file, an HDF5-compatible file format. We will also save the latitude and longitude of the points on the sphere.

``` @example DEC
sim = eval(gensim(budyko_sellers_halfar_water))
fₘ = sim(s, generate)

tₑ = 100.0

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

save("ice.jld2",
  Dict("lat" => map(x -> x[1], p_sph), "lon" => map(x -> x[2], p_sph), "ice" => soln(tₑ).h))

(extrema(soln(0.0).h), extrema(soln(tₑ).h))
```

## Visualize 

Let's visualize the initial conditions for ice height and the ice height after 100 years.

### Initial ice height

``` @example DEC
f = Figure()
ax = LScene(f[1,1], scenekw=(lights=[],))
update_cam!(ax.scene, Vec3f(0,0,0.8), Vec3f(0,0,0), Vec3f(0, 1, 1))
msh = mesh!(ax, s_plots, color=soln.u[begin].h, colormap=Reverse(:redsblues))
Colorbar(f[1,2], msh)
f
```

### Ice height after 100 years

``` @example DEC
f = Figure()
ax = LScene(f[1,1], scenekw=(lights=[],))
update_cam!(ax.scene, Vec3f(0,0,0.8), Vec3f(0,0,0), Vec3f(0, 1, 1))
msh = mesh!(ax, s_plots, color=soln.u[end].h, colorrange=extrema(soln.u[begin].h), colormap=Reverse(:redsblues))
Colorbar(f[1,2], msh)
f
```

```@example DEC
run(`rm $ice_thickness_file`) # hide
```

