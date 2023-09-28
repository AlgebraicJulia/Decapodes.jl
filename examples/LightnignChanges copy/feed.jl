"""
The FEED (Fluid Energetic Electron Deposition) model computes electron
deposition profiles given incoming electron flux. We can model these dynamics on
a vertical 1-manifold embedded in the 3D ionosphere.

One could say that the question an electron deposition model seeks to answer is
"Where are the electrons deposited?" We can answer this with two key related
quantities.

First, we can provide the "penetration depth", R. This is the altitude at which
enough energy has dissipated that we consider an electron's effects on its
neighbors to be neglible. This is like measuring where a round entering
ballistic gel stops.

Second, we can provide a measure of "amount of electrons deposited" by keeping
track of how electrons are freed upon colliding with species in the ionosphere. Electrons of higher energy hit more neutral
species, freeing many electrons, but are few in number. Electrons of low energy
collide with fewer molecules, but are more common. Electrons colliding with other species to produce more electrons is called "secondary ionization."

The "primary" ionization is the result of electrons precipitating into our
atmosphere. The WIPP (Whistler-Induced Particle Precipitation) model gives the
incoming electron flux at 200 km. As a pre-processing step, FEED takes the
single value of electron flux WIPP gives at 200km, and "integrates [a
differential equation derived by Kim] progressively downwards" to give values
for ϕₑ across the domain.

Kim's differential equation describes how the ionosphere resists electron flux.
The term that informs resistance is known as "stopping force", and is given by
Berger et al. in a quantized relation between energy and stopping force, which
we can interpolate.

Given electron flux, ϕₑ, we would like to compute the trail of energy that these
electrons leave. Coulomb scattering renders it infeasible to model electron
collisions directly. So, scientists like Gruen and Rees have fit statistical
power laws to predict penetration depth. Maeda (1964) succinctly describes such
models as "semi-empirical".

The engine of FEED that computes energy across the domain is a differential
equation which describes ionization rate at discrete altitudes. This
differential equation, given by Rees, relies on a power law describing
penetration depth, given by Gledhill. Thus, our first answer to our question, R,
informs the second answer, E.

FEED is somewhere in the middle of the spectrum of empirical vs. statistical
models.

The components of the FEED model can be organized like so:
- A differential equation (Eq. 2-16 Kim 2020), "Kim", describing ∂ϕₑ/∂z is integrated to provide values for electron flux, ϕₑ, across the domain.
  - This differential equation depends on an estimation of electron stopping force, F, given by Berger et al. 1984. Call this component "Berger".
- A differential equation (Eq. 1 Rees 1963) computes ionization rate. Call this component "Rees".
  - This differential equation depends on a power law describing penetration depth (Eqs. 5 & 6 Gledhill 1973). Call this component "Gledhill".

The physical quantities to be tracked here are thus:
- R
  - The penetration depth, usually measured in pressure altitude.
- E
  - Energy. In practice, we want to measure how many electrons are produced at pre-quantized energy levels. Typically at 10 keV, 50 keV, 100 keV, and so on.


So why use Decapodes for FEED?
Decapodes excels at switching out component physics, since composition patterns are treated as first-class objects. This is especially relevant here, since the power law to use for R has been hotly debated by Gruen, Rees, Gledhill, Spencer, Maeda, and others. See Gledhill's introduction for a short history of this debate.
Decapodes decouples the compiled computation code from the model representation. This complements our first advantage, in that once we explicitly represent the physics of our new component, we do not need to write its solver by hand.
The way that models are written in a Decapode more closely matches the way that an equation is written on a blackboard, or how it appears in a paper. When the model and the solver are mixed in some imperative way, the model is not readily visible. Given the code alone, it is usually not clear what the model is at all.
The implementation of a model in code influences the way that scientists represent their work in the literature. We note that in the Kim dissertation, there is a "hidden" differential equation. There is a paragraph describing how certain quantities are multiplied by a discrete timestep and later summed. What is being described is a explicit timestepping solver of a differential equation. Thus, the compiled code is translated into English without ever appearing as a differential equation!

##############
# References #
##############

# * Denotes that an equation from one of these papers is currently directly
# modeled.

# *
# Kim, D. (2020), Analysis of conjugate lightning-induced electron precipitation
#   events using the VLF remote sensing method, Ph.D. dissertation, University
#   of Florida, Gainesville, Florida. 

# *
# Rees, M. H., “Auroral ionization and excitation by incident energetic
#   electrons”, Planetary and Space Science, vol. 11, no. 10, pp. 1209–1218,
#   1963. doi:10.1016/0032-0633(63)90252-6.

# *
# Gledhill, J. A., “The range-energy relation for 0.1-600 keV electrons,”
#   Journal of Physics A: Mathematical, Nuclear and General, vol. 6, no. 9, pp.
#   1420–1428, 1973. doi:10.1088/0305-4470/6/9/017 

# Maeda, K., “Diffusion of low energy auroral electrons in the atmosphere,”
#   Journal of Atmospheric and Terrestrial Physics, vol. 27, no. 2, pp. 259–275,
#   1965. doi:10.1016/0021-9169(65)90121-2 

# Berger M. J., Seltzer S. M. "Stopping powers and ranges of electrons and
#   positrons , 2nd Ed. U.S. Department of Commerce, National Bureau of
#   Standards. 1983
"""

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
#using GLMakie
using CairoMakie
using Interpolations
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
s = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s, 20001, point=Point2D.(0,0:1e-2:200))
add_edges!(s, 1:(nv(s)-1), 2:nv(s))
sd = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s)
subdivide_duals!(sd, Circumcenter())
(nv(sd),ne(sd))
end # Load the Mesh

#############
# Operators #
#############
begin # Operators

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :integrate => (x,y) -> begin
      fill(x, nv(sd))
    end
    :log10 => x -> begin
      log10.(x)
    end
    :.* => (x,y) -> x .* y
    :./ => (x,y) -> x ./ y
    :.+ => (x,y) -> x .+ y
    :.- => (x,y) -> x .- y
    :.^ => (x,y) -> x .^ y
    :^ => (x,y) -> x ^ y
    # TODO: Will this shorthand cause problems when we are trying to subtract?
    :- => x -> -1 * x
    #:- => (x,y) -> x .- y
    :.> => (x,y) -> 1 * (x .> y)
    :.< => (x,y) -> 1 * (x .< y)
    :.≤ => (x,y) -> 1 * (x .≤ y)
    :.≥ => (x,y) -> 1 * (x .≥ y)

    #_ => error("Unmatched operator $my_symbol")
    #_ => default_dec_generate(sd, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end
end # Operators

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

# TODO: Let E be a constant for now.
# To get physical synthetic data, you can model E = a*sin(c*t) + b. s.t. b>a.

constants_and_parameters = (
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
# Interpolation #
#################

begin
using Interpolations
using CairoMakie
# TODO: Instead of interpolating, can we just reproduce the table from Eq.s 1 &
# 2 from Rees?
z = [2.34e-1, 1.37e-1, 7.76e-2, 2.65e-7]
h = [60, 64, 68, 300]
h_new = 60:.2:200
interp_linear = linear_interpolation(h, z)
z_new = map(p -> interp_linear(p), h_new)
end

####################
# Component Models #
####################
begin # Component Model

# Gruen 1957 uses this simple power law for R, for values in [5, 54] keV.
Gruen = @decapode begin
  R::Form0

  R == 4.57e-6 * E ^ 1.75
end
# Rees 1963 argued that it is still rather physical even at 0.4 keV and 300 keV.
# Rees later revised this, as Gledhill 1973 notes.

# Eqs. 2-18 & 2-19 from Kim 2020
# Originally Eqs. 5 & 6 from Gledhill 1973.
Gledhill = @decapode begin
  (R,Z)::Form0
  E::Parameter
  n_M::Constant

  R == (E .≤ 100) .* (10^(-5.133 + 1.358*log10(E)+0.215*(log10(E))^2-0.043*(log10(E))^3 )) +
       (E .> 100) .* (10^( -6.193+2.585*log10(E)-0.22*(log10(E))^2))
  # TODO: Why quantize? Can't we just interpolate?
  # TODO: Further, is there even a need for interpolation if we just compute the
  # "ionization production" using Eq.s 1 & 2 from Rees?
  n_M_R == find_closest_n_M_to_R(n_M, R)
end

# TODO: Integrate ϕ given ϕ at a single point in equation 2-16.
Kim = @decapode begin
  ϕₑ::Form0
  (E, ϕ₀)::Parameter

  # TODO: Represent more of this integration explicitly.
  # TODO: Just use numerical integration given Eq. 2-16.
  ϕₑ == integrate(ϕ₀, E)
end

Rees = @decapode begin
  (ΔEᵢₒₙ, ρ, n_M)::Constant
  (ionization, R, Z, ϕₑ)::Form0
  E::Parameter # Form0

  r₀ == R / ρ

  # TODO: Energy dissipation function, λ comes from a reference table, from Spencer 1959 pp.38.
  # q is rate of ionization.

  # Eq. 2-17 from Kim 2020, and Eq. 1 from Rees 1963.
  q == ∂ₜ(ionization)
  q == ((E / r₀) / ΔEᵢₒₙ * λ(Z/R) * (n_M)/(n_M_R)) * ϕₑ
end
#to_graphviz(Rees)

end # Component Model

#####################
# Model Composition #
#####################
begin # Model Composition
draw_composition_diagram(dgm) = to_graphviz(dgm, junction_labels=:variable, box_labels=:name, prog="circo")

# Compose the dynamics for a single "species".
compose_precipitation = @relation () begin
  Rees(E, ϕₑ, R, Z, n_M, n_M_R)
  Kim(E, ϕₑ)
  Gledhill(R, E, Z, n_M, n_M_R)
end
draw_composition_diagram(compose_precipitation)

precipitation_cospan = oapply(compose_precipitation, [
   Open(Rees,     [:E, :ϕₑ, :R, :Z, :n_M, :n_M_R]),
   Open(Kim,      [:E, :ϕₑ]),
   Open(Gledhill, [:R, :E, :Z, :n_M, :n_M_R])])
Precipitation = apex(precipitation_cospan)
infer_types!(Precipitation, op1_inf_rules_1D, op2_inf_rules_1D)
to_graphviz(Precipitation)

# Compose the dynamics for a single "species".
# Note that I construct this composition pattern "imperatively" here.
# There is an example using a "declarative" approach at the bottom of this file.
function quantize_multi_species(quants = [50, 100, 200, 300, 500])
  dgm = RelationDiagram(0)
  add_parts!(dgm, :Box, length(quants); name=map(quants) do q
      Symbol("E_"* replace(string(q), r"\." => "_"))
    end)
  add_parts!(dgm, :Junction, 2;
    variable=[:Z, :n_M])
  add_parts!(dgm, :Port, length(quants)*2;
    box=repeat(1:length(quants), inner=2),
    junction=repeat(1:2, length(quants)))
  dgm, apex(oapply(dgm, 
    fill(Open(Precipitation, [:Z, :n_M]), length(quants))))
end
compose_MSP, MSP = quantize_multi_species()
draw_composition_diagram(compose_MSP)

MSP = expand_operators(MSP)
infer_types!(MSP, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(MSP, op1_res_rules_1D, op2_res_rules_1D)

to_graphviz(MSP)
end # Model Composition

#########################
# Simulation Generation #
#########################
begin # Simulation Generation

open("./generated_precipitation_sim.jl", "w") do file
  write(file, string(gensim(MSP)))
end
sim = include(eval, "../../generated_precipitation_sim.jl")
fₘ = sim(sd, generate)

end # Simulation Generation

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
   #map(x -> Symbol(x == :e ? :ρ_e : "Chemistry_ρ_" * string(x)), species)...,
   map(x -> Symbol(x == :e ? :ρ_e : "Chemistry_ρ_" * string(x)), collect(keys(densities)))...,
   #map(x -> Symbol("Chemistry_" * string(x)), rateCoef_names)...])
   map(x -> Symbol("Chemistry_" * string(x)), collect(keys(rateCoefs)))...])

tₑ = 2.0 # [s]
#dt = 20e-3 [s]

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5(), progress=true, progress_steps=1)

@save "feed.jld2" soln

############
# Old Code #
############

################################
# Non-Hierarchical Composition #
################################
#compose_multi_species = @relation () begin
#  Rees50(E_50, ϕₑ_50, R_50, Z, n_M, n_M_R_50)
#  Kim50(E_50, ϕₑ_50)
#  Gledhill50(R_50, E_50, Z, n_M, n_M_R_50)
#
#  Rees100(E_100, ϕₑ_100, R_100, Z, n_M, n_M_R_100)
#  Kim100(E_100, ϕₑ_100)
#  Gledhill100(R_100, E_100, Z, n_M, n_M_R_100)
#
#  Rees200(E_200, ϕₑ_200, R_200, Z, n_M, n_M_R_200)
#  Kim200(E_200, ϕₑ_200)
#  Gledhill200(R_200, E_200, Z, n_M, n_M_R_200)
#
#  Rees300(E_300, ϕₑ_300, R_300, Z, n_M, n_M_R_300)
#  Kim300(E_300, ϕₑ_300)
#  Gledhill300(R_300, E_300, Z, n_M, n_M_R_300)
#
#  Rees500(E_500, ϕₑ_500, R_500, Z, n_M, n_M_R_500)
#  Kim500(E_500, ϕₑ_500)
#  Gledhill500(R_500, E_500, Z, n_M, n_M_R_500)
#end
#to_graphviz(compose_multi_species, junction_labels=:variable, box_labels=:name, prog="circo")
#
#multi_species_precipitation_cospan = oapply(compose_multi_species, [
#   Open(Rees,     [:E, :ϕₑ, :R, :Z, :n_M, :n_M_R]),
#   Open(Kim,      [:E, :ϕₑ]),
#   Open(Gledhill, [:R, :E, :Z, :n_M, :n_M_R]),
#
#   Open(Rees,     [:E, :ϕₑ, :R, :Z, :n_M, :n_M_R]),
#   Open(Kim,      [:E, :ϕₑ]),
#   Open(Gledhill, [:R, :E, :Z, :n_M, :n_M_R]),
#
#   Open(Rees,     [:E, :ϕₑ, :R, :Z, :n_M, :n_M_R]),
#   Open(Kim,      [:E, :ϕₑ]),
#   Open(Gledhill, [:R, :E, :Z, :n_M, :n_M_R]),
#
#   Open(Rees,     [:E, :ϕₑ, :R, :Z, :n_M, :n_M_R]),
#   Open(Kim,      [:E, :ϕₑ]),
#   Open(Gledhill, [:R, :E, :Z, :n_M, :n_M_R]),
#
#   Open(Rees,     [:E, :ϕₑ, :R, :Z, :n_M, :n_M_R]),
#   Open(Kim,      [:E, :ϕₑ]),
#   Open(Gledhill, [:R, :E, :Z, :n_M, :n_M_R])
#   ])

############################
# Hierarchical Composition #
############################
#compose_multi_species = @relation () begin
#  E_50( Z, n_M)
#  E_100(Z, n_M)
#  E_200(Z, n_M)
#  E_300(Z, n_M)
#  E_500(Z, n_M)
#end
#draw_composition_diagram(compose_multi_species)
#
#multi_species_precipitation_cospan = oapply(compose_multi_species, [
#   Open(Precipitation,     [:Z, :n_M]),
#   Open(Precipitation,     [:Z, :n_M]),
#   Open(Precipitation,     [:Z, :n_M]),
#   Open(Precipitation,     [:Z, :n_M]),
#   Open(Precipitation,     [:Z, :n_M])])
#MSP = apex(multi_species_precipitation_cospan)
