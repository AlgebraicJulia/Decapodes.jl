##############
# References #
##############

# F. Heidler, J. M. Cvetic and B. V. Stanic, "Calculation of lightning current
#   parameters," in IEEE Transactions on Power Delivery, vol. 14, no. 2, pp.
#   399-404, April 1999, doi: 10.1109/61.754080.

# Kotovsky, D. A. (2016), Response of the nighttime upper mesosphere to electric
#   field changes produced by lightning discharges, Ph.D. dissertation,
#   University of Florida, Gainesville, Florida. 

# V. P. Pasko, U. S. Inan, T. F. Bell, and Y. N. Taranenko, “Sprites produced by
#   quasi-electrostatic heating and ionization in the lower ionosphere,”
#   Journal of Geophysical Research: Space Physics, vol. 102, no. A3,
#   pp. 4529–4561, 1997, doi: 10.1029/96JA03528.

# G. Veronis, V. P. Pasko, and U. S. Inan, “Characteristics of mesospheric
#   optical emissions produced by lightning discharges,” Journal of Geophysical
#   Research: Space Physics, vol. 104, no. A6, pp. 12645–12656, 1999,
#   doi: 10.1029/1999JA900129.

################
# Dependencies #
################
using Pkg
Pkg.add("Interpolations")
Pkg.add("Roots")
Pkg.add("MAT")
Pkg.add("CairoMakie")
Pkg.add("TerminalLoggers")
Pkg.resolve()
begin # Dependencies
# AlgebraicJulia
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes

# External Dependencies
using Base.MathConstants: e, π
using Distributions
using CairoMakie
using Interpolations
using JLD2
using LinearAlgebra
using Logging
using MAT
using MLStyle
using MultiScaleArrays
using OrdinaryDiffEq
using Roots
using SparseArrays
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}
end # Dependencies

#################
# Load the Mesh #
#################
begin # Load the Mesh
# We assume cylindrical symmetry.
MAX_r = 400.0e3 # [m]
MAX_Z = 100.0e3 # [m]

include("../grid_meshes.jl")
s = triangulated_grid(400.0e3, 100.0e3, 500e1, 500e1, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())
nv(s), ne(s), ntriangles(s)
wireframe(s, linewidth=1/2)
wireframe(sd, linewidth=1/2)

end # Load the Mesh

#############
# Operators #
#############
begin # Operators

include("./operators.jl")

function generate(sd, my_symbol; hodge=GeometricHodge())
  star0_mat = ⋆(0,sd,hodge=hodge)
  star1_mat = ⋆(1, sd, hodge=hodge)
  star2_mat = ⋆(2, sd, hodge=hodge)
  invstar0_mat = inv_hodge_star(0, sd, hodge=hodge)
  invstar1_mat = inv_hodge_star(1, sd, hodge=hodge)
  d1_mat = d(1,sd)
  duald0_mat = dual_derivative(0,sd)
  avg_mat = avg₀₁(sd)
  below_60_mask = map(x -> x[2] < 60*1000, sd[:point])
  sharp_mat = ♯_mat(sd, DiscreteExteriorCalculus.AltPPSharp())
  op = @match my_symbol begin
    :⋆ => x -> begin
      error("Uninferred hodge star")
    end
    :d => x -> begin
       error("Uninferred exterior derivative")
    end
    :⋆₀ => x -> begin
      # We create a local mat variable as a hint for Julia to cache it.
      star0_mat*x
    end
    :⋆₁ => x -> begin
      star1_mat*x
    end
    :⋆₂ => x -> begin
      star2_mat*x
    end
    :⋆₀⁻¹ => x -> begin 
      invstar0_mat*x
    end
    :⋆₁⁻¹ => x -> begin
      invstar1_mat*x
    end
    :d₁ => x -> begin
      d1_mat*x
    end
    :d̃₀ => x -> begin
      duald0_mat*x
    end
    :avg₀₁ => x -> begin
      avg_mat*x
    end
    :exp => x -> begin
      exp.(x)
    end
    :sqrt => x -> begin
      sqrt.(x)
    end
    :∧ᵖᵈ => (x,y) -> begin
      pd_wedge(Val{(1,1)}, sd, x, y)
    end
    :log10 => x -> begin
      log10.(x)
    end
    :exp10 => x -> begin
      exp10.(x)
    end
    :maskalt => x -> begin
      x[below_60_mask] .= 0.0
      x
    end
    :♯ => x -> begin
      sharp_mat * x
    end
    :mag => x -> norm.(x)
    :.* => (x,y) -> x .* y
    :./ => (x,y) -> x ./ y
    :.+ => (x,y) -> x .+ y
    :.- => (x,y) -> x .- y
    :.^ => (x,y) -> x .^ y
    :^ => (x,y) -> x ^ y
    :- => x -> -1 * x
    #:- => (x,y) -> x .- y
    :.> => (x,y) -> 1 * (x .> y)
    :.< => (x,y) -> 1 * (x .< y)
    :.≤ => (x,y) -> 1 * (x .≤ y)
    :.≥ => (x,y) -> 1 * (x .≥ y)
    :invert_mask => x -> (!).(x)
    :id => x -> x

    #_ => error("Unmatched operator $my_symbol")
    _ => default_dec_generate(sd, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end
end # Operators

######################
# Heidler Parameters #
######################
begin # Heidler Parameters

T1 = 50.0         # [us] See Kotovsky Table 5-1 column 1
T2 = 1000.0       # [us]
τ₁ = T1 * 1e-6    # [s]
τ₂ = T2 * 1e-6    # [s]
n = 10.0          # Heidler normalizing constant

# TODO: Should this distribution should reach its max at 50 microseconds?
# This is equation 5-14 from Kotovsky, or 8 from Heidler.
# Heidler uses an intermediate variable kₛ = t/τ₁.
HeidlerForη = t -> (t / τ₁)^n / (1 + ((t /τ₁)^n)) * exp(-1.0 * t / τ₁)
ddtHeidlerForη = t -> (exp(-t / τ₁) * (t/τ₁)^n * ((t*(-(t/τ₁)^n - 1)) + n*τ₂)) /
  (t * τ₂ * ((t/τ₁)^n + 1)^2)
time_for_k = find_zero(ddtHeidlerForη, (1.0e-6, 100.0e-6)) # 8.0878657...e-5
#plot(0.0:1.0e-8:100.0e-6, HeidlerForη)
#log_Heidler = map(0.0:1.0e-5:1.5e-3) do x
#  log(HeidlerForη(x)*1e6)
#end
#log_Heidler[1] = log_Heidler[2] # Set this manually to something not -Inf
#plot(0.0:1.0e-5:1.5e-3, log_Heidler)
#plot(0.0:1.0e-8:100.0e-6, ddtHeidlerForη)
η = HeidlerForη(time_for_k) # Normalizing constant 0.196775...

end # Heidler Parameters

#############
# Constants #
#############
begin # Constants
#temp = (return_time ./ rise_time) .^ 10 # Relative time difference.
I₀ = 250.0e3      # [kA] See Kotovsky pp. 107
ε₀ = 8.854e-12    # [F m⁻¹] Permittivity of free space
μ₀ = 4*π*1e-7     # [H m⁻¹] Permeability of free space
Z₀ = sqrt(μ₀/ε₀)  # [Ω] Impedance of free space
c = sqrt(1/ε₀/μ₀) # [m s⁻¹] Speed of light
v = 2/3 * c       # [m s⁻¹] Approximate speed of propagation of lightning strike through medium
qₑ = 1.602e-19    # [coul] Charge of electron
mₑ = 9.109e-31    # [kg] Mass of electron
kB = 8.617e-5     # [eV K⁻¹] Boltzmann constant
# Note that these ρ₀ and z₀ are chosen for numerical stability.
# For a rectangular grid discretization, good values are 3Δρ, 3Δz.
z₀ = 3.0e3        # [km] Heidler normalizing constant
ρ₀ = 3.0e3        # [km] Heidler normalizing constant
a = 10.0 * 1e3    # [m] Starting altitude of vertical decay for Heidler
# TODO: If the axes are swapped, then these accessors must be as well.
Z = EForm(avg₀₁(sd) * VForm(map(p -> p[2], sd[:point]))) # [m] Altitude
ρ = EForm(avg₀₁(sd) * VForm(map(p -> p[1], sd[:point]))) # [m] Distance from centre of strike

constants_and_parameters = (
  #mₑ = mₑ,
  #μ₀ = μ₀,
  #Z₀ = Z₀,
  #e = e,
  Chemistry_kB = kB,
  Veronis_qₑ = qₑ,
  Veronis_c = c,
  Veronis_ε₀ = ε₀,
  Heidler_τ₁ = τ₁,
  Heidler_τ₂ = τ₂,
  Heidler_I₀ = I₀,
  Heidler_v = v,
  Heidler_n = n,
  Heidler_a = a,
  Heidler_η = η,
  Heidler_t = t -> t,
  Heidler_z₀ = z₀,
  Heidler_π = π,
  Heidler_ρ₀ = ρ₀)

end # Constants

#################
# Heidler Model #
#################
begin # Heidler Model

# This defines J.
# The model takes in rise time, fall time, and peak current of the lightning
# strike, (i₀).
# Note that this model does not have any tangent variables.

# See Veronis et al. § 2.

Heidler = @decapode begin
  (J, J₀, Z, ρ, tret)::Form1{X}
  (τ₁, τ₂, I₀, v, n, a, η, z₀, π, ρ₀)::Constant{X}
  (t)::Parameter{X}

  tret == t .- Z ./ v
  temp == (tret./τ₁).^n
  # See Veronis et al. § 2.1
  J₀ == (1.0 / (π * ρ₀^2)) * I₀/η * temp ./ (1 .+ temp) .* exp(-1 * tret ./ τ₂) .* (tret .> 0)

  # See Kotovsky Eq. 5-11a and 5-11b.
  J == (Z .≤ a) .* J₀ .* exp(-1.0 .* ρ .^ 2 / ρ₀ .^ 2) +
       (Z .> a) .* J₀ .* exp((-1.0 * ρ .^ 2 / ρ₀^2) - ((Z .- a) .^ 2 / z₀^2))
end
to_graphviz(Heidler)

end # Heidler Model

###########
# Veronis #
###########
begin # Veronis
# Assumptions that allow for cylindrical "pseudo-3D":
# - Lightning can be restricted to plane.
# - The magnetic field is orthogonal to this plane.

# More assumptions:
# - Tn is the neutral temperature. We assume that this never changes due to its
#   magnitude, and the timescale of investigation being small. As a consequence,
#   it acts as something like an "infinite source."

Veronis = @decapode begin
  # TODO: Double check units on these (esp. w.r.t. σ being integrated along
  # edges.)
  B::DualForm0{X}
  E::Form1{X}
  θ::Form0{X}
  J::Form1{X}
  (ρ_e, ρ_gas, Tn)::Form0{X}
  σ::Form1{X}
  (qₑ,c,ε₀)::Constant{X}

  # See Veronis et al. Equations 1 & 2
  Ė == -1 * (J - σ .* E)./ε₀ + ((c^2) .* ⋆(d(B)))
  Ė == ∂ₜ(E)

  # See Veronis et al. Equation 3
  Ḃ == ⋆(d(E))
  Ḃ == ∂ₜ(B)

  # See Kotovsky pp. 91: theta = E/n
  θ == (1e21/1e6) .* mag(♯(E)) ./ ρ_gas

  Eq5_2a_mask == .>(θ, (0.0603 * sqrt(200 ./ Tn)))
  Eq5_2b_mask == .≤(θ, (0.0603 * sqrt(200 ./ Tn)))
  # See Kotovsky pp. 92 5-2a
  Eq5_2a == qₑ * ρ_e ./ ρ_gas *
    (10 .^ ( 50.97 + (3.026 * log10( 1e-21*θ )) + (8.4733e-2 * (log10( 1e-21*θ ) .^ 2))))
  # See Kotovsky pp. 92 5-2b
  Eq5_2b == qₑ * 3.656e25 * ρ_e ./ ρ_gas .* sqrt(200.0 ./ Tn)

  σ == avg₀₁((Eq5_2a_mask .*  Eq5_2a) + (Eq5_2b_mask .*  Eq5_2b))
end
#to_graphviz(Veronis)
#to_graphviz(resolve_overloads!(infer_types!(Veronis)))
end # Veronis

################################################
# Initialize Densities and Neutral Temperature #
################################################
begin # Format Atmosphere

include("formatAtmosphere.jl")
species, densities, Tn, rateCoef_names, rateCoefs = formatAtmosphere("./examples/LightnignChanges copy/chi180_O_dyn.mat", sd)
mesh(s, color=densities[:O])

end # Format Atmosphere

#############################
# Initialize E, B, σ, and θ #
#############################
begin # Initialize Electromagnetism

Ef = VForm(zeros(ne(s)))
B = DualForm{0}(zeros(ntriangles(s)))

# Initialize reduced electric field, θ [V m m]
θ = nothing
E₀ = nothing
# Linear model of the reduced electric field
#if option == 3
  E₀ = Ef
  #θ = 1e21/1e6 * E ./ avg₀₁(sd, densities[:gas]) # Danny Dissert p 91 theta = E/n
  θ = EForm(zeros(ne(s)))
#end

# Initialize conductivity
#σ = zeros(nv(s))
for p in vertices(s)
  # Note: This is assuming that points in s are in Cartesian coordinates.
  ionosphere_altitude = 60.0e3
  #if s[p, :point][3] < ionosphere_altitude
  if s[p, :point][2] < ionosphere_altitude
    #σ[p] = 0.0
  else
    # See Kotovsky pp.92 eq 5-2b
    # See also Pasko et al.
    #σ[p] = qₑ * 3.656e25 * densities[:e][p] / densities[:gas][p] * sqrt(200.0 / Tn[p])
  end
end
#mesh(s, color=σ)

end # Initialize Electromagnetism

#####################
# Model Composition #
#####################
begin # Model Composition

# TODO: Rename this file to chemistry.jl
include("rateCoefficients_dynam.jl")

compose_lightning = @relation () begin
  Heidler(J)
  Veronis(θ, J, ρ_e, ρ_gas, Tn)
  Chemistry(Tn, θ, ρ_e, ρ_gas)
end
#to_graphviz(compose_lightning, junction_labels=:variable, box_labels=:name, prog="circo")

lighting_cospan = oapply(compose_lightning,
  [Open(Heidler, [:J]),
  Open(Veronis, [:θ, :J, :ρ_e, :ρ_gas, :Tn]),
  Open(Chemistry, [:Tn, :θ, :ρ_e, :ρ_gas])])

lightning = apex(lighting_cospan)
# Warning: This diagram is large.
#to_graphviz(lightning)

lightning = ∘(resolve_overloads!, infer_types!, expand_operators)(lightning)

end # Model Composition

#########################
# Simulation Generation #
#########################
begin # Simulation Generation

#open("./generated_lightning_sim.jl", "w") do file
#  write(file, string(gensim(lightning)))
#end
sim = include("../../generated_lightning_sim.jl")
fₘ = sim(sd, generate)

end

###########
# Solving #
###########

u₀ = construct(PhysicsState,
  map(x -> VectorForm(x.data),
    [Z,
     ρ,
     B,
     Ef,
     Tn,
     values(densities)...,
     values(rateCoefs)...]),

  Float64[],

  [:Heidler_Z,
   :Heidler_ρ,
   :Veronis_B,
   :Veronis_E,
   :Tn,
   map(x -> Symbol(x == :e ? :ρ_e : "Chemistry_ρ_" * string(x)), collect(keys(densities)))...,
   map(x -> Symbol("Chemistry_" * string(x)), collect(keys(rateCoefs)))...])

# This is a typical value for a model run:
#tₑ = 200e-6 # [s]

# This is how long it takes light to travel 400 km in a vacuum.
#tₑ = 1.334e-3 # [s] 

tₑ = 200e-6 * .67
prob = ODEProblem{true}(fₘ, u₀, (0, tₑ), constants_and_parameters)
@time soln = solve(prob, Tsit5(), progress=true, progress_steps=1, save_everystep=false)

# Warning: This overwrites a file with this name if there is one.
@save "lightning.jld2" soln

# There is no need to load the solution back in if you haven't killed the REPL.
@load "lightning.jld2" soln 

############
# Plotting #
############

# Compute the divergence. i.e. ⋆d⋆ in the Discrete Exterior Calculus
div(form1) = inv_hodge_star(0,sd)*dual_derivative(1,sd)*hodge_star(1,sd)*form1

# Plot the divergence of the E field. i.e. ⋆d⋆(E)
function plot_div_E(t, colorrange=extrema(findnode(soln(t), :Veronis_E)))
  mesh(s, color=div(findnode(soln(t), :Veronis_E)); colorrange=colorrange)
end

plot_div_E(tₑ)
plot_div_E(tₑ, (-1e-32, 1e-32))

# Plot the difference in electron density at each point.
function plot_diff_e(t)
  mesh(s, color=findnode(soln(t), :ρ_e) - findnode(soln(0), :ρ_e))
end

plot_diff_e(tₑ)

# Plot the density of e over the initial density of e.
function plot_px_e(t)
  mesh(s, color=findnode(soln(tₑ), :ρ_e) / findnode(soln(0), :ρ_e))
end

plot_px_e(tₑ)

# Plot the difference in O density at each point.
function plot_diff_O(t)
  mesh(s, color=findnode(soln(tₑ), :Chemistry_ρ_O) - findnode(soln(0), :Chemistry_ρ_O))
end

plot_diff_O(tₑ)

# Plot the density of O over the initial density of O.
function plot_px_O(t)
  mesh(s, color=findnode(soln(tₑ), :Chemistry_ρ_O) / findnode(soln(0), :Chemistry_ρ_O))
end

plot_px_O(tₑ)

# Plot the E field itself.
function plot_E(t)
  sharp_mat = ♯_mat(sd, AltPPSharp())
  CairoMakie.arrows(point(s), sharp_mat*findnode(soln(t), :Veronis_E), lengthscale=1e3)
end

plot_E(tₑ)

# Plot the magnitude of the E field.
function plot_mag_E(t)
  sharp_mat = ♯_mat(sd, AltPPSharp())
  mesh(s, color=norm.(sharp_mat*findnode(soln(t), :Veronis_E)))
end

plot_mag_E(tₑ)

# Plot the log of the absolute value of the B field.
function plot_log_abs_B(t)
  scatter(dual_point(sd, triangle_center(sd, triangles(sd))),
    color=(log ∘ abs).(findnode(soln(t), :Veronis_B)))
end

plot_log_abs_B(tₑ)

# Plot the B field itself.
function plot_B(t)
  scatter(dual_point(sd, triangle_center(sd, triangles(sd))),
    color=findnode(soln(t), :Veronis_B))
end
plot_B(tₑ)

# Find the largest nonzero value of the given 0Form.
function find_max_x_nonzero_form0(form)
  max_x = 0
  for (i,p) in enumerate(point(sd))
    if form[i] != 0
      max_x = p[1]
    end
  end
  max_x
end
find_max_x_nonzero_form0(findnode(soln(tₑ), :ρ_e))

# Find the largest nonzero value of the given 1Form.
function find_max_x_nonzero_form1(form)
  max_x = 0
  for i in edges(sd)
    p = point(sd, s[i, :∂v0])
    if form[i] != 0
      max_x = p[1]
    end
  end
  max_x
end
find_max_x_nonzero_form1(findnode(soln(tₑ), :Veronis_E))

# Create a gif of magnitude of E.
# This requires that you did not set save_everystep to false when you solved.
function make_gif_E()
  frames = 100
  ♯_m = ♯_mat(sd, DiscreteExteriorCalculus.AltPPSharp())
  fig, ax, ob = mesh(s, color=norm.(♯_m * findnode(soln(0), :Veronis_E))), colormap=:jet, colorrange=extrema(norm.(♯_m * findnode(soln(tₑ), :Veronis_E)))
  Colorbar(fig[1,2], ob)
  record(fig, "lightning.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    ob.color = norm.(♯_m * findnode(soln(t), :Veronis_E))
  end
end

make_gif_E()
