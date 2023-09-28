# Use @inbounds.
# Try sum([...]), instead of (.+)(...).
# Try vcat(...)*ones(nv(s))
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
using Catlab, Catlab.Graphics
using CombinatorialSpaces, CombinatorialSpaces.ExteriorCalculus
using Decapodes

# External Dependencies
using Distributions, CairoMakie, Interpolations, JLD2, LinearAlgebra, Logging
using MAT, MLStyle, MultiScaleArrays, OrdinaryDiffEq, Roots, SparseArrays
using Base.MathConstants: e, π
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
MAX_r = 400.0e3 # [m]
MAX_Z = 100.0e3 # [m]
include("./examples/grid_meshes.jl")

# "Our numerical simulation box extends from z=0 km to z=100 km and from r=0 km
# to r=400 km, with grid spacings of Δz=0.5 km and Δr=2 km."
# ~ Veronis § 2.1, pp. 12646
# Veronis Grid Resolution, Triangulated:
# TODO: See effect of declaring the mesh `const`.
s = triangulated_grid(400.0e3, 100.0e3, 0.5e3, 2e3, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())
nv(s), ne(s), ntriangles(s)
end # Load the Mesh

#############
# Operators #
#############
begin # Operators

# Definitions of avg₀₂, and ∧₁₁′:
include("./examples/LightnignChanges copy/operators.jl")

# These are the matrices of the differential operators.
hodge=DiagonalHodge();
star0_mat = ⋆(0,sd,hodge=hodge);
star1_mat = ⋆(1, sd, hodge=hodge);
star2_mat = ⋆(2, sd, hodge=hodge);
invstar0_mat = inv_hodge_star(0, sd, hodge=hodge);
invstar1_mat = inv_hodge_star(1, sd, hodge=hodge);
d1_mat = d(1,sd);
duald0_mat = dual_derivative(0,sd);
avg_mat = avg₀₁(sd);
below_60_mask = map(x -> x[2] < 60*1000, sd[:point]);
sharp_mat = ♯_mat(sd, DiscreteExteriorCalculus.AltPPSharp());

function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :⋆ => x -> begin
      error("Uninferred hodge star")
    end
    :d => x -> begin
       error("Uninferred exterior derivative")
    end
    :⋆₀ => x -> begin
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
    #:maskneg => x -> begin
    #  x[x .< 0.0] .= 0.0
    #  x
    #end
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
    :.> => (x,y) -> 1 * (x .> y)
    :.< => (x,y) -> 1 * (x .< y)
    :.≤ => (x,y) -> 1 * (x .≤ y)
    :.≥ => (x,y) -> 1 * (x .≥ y)
    :invert_mask => x -> (!).(x)
    :id => x -> x

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
η = HeidlerForη(time_for_k) # Normalizing constant 0.196775...

end # Heidler Parameters

#############
# Constants #
#############
begin # Constants
#I₀ = 250.0e3      # [kA] See Kotovsky pp. 107
I₀ = 250.0e3      # [A] See Kotovsky pp. 107
ε₀ = 8.854e-12    # [F m⁻¹] Permittivity of free space
μ₀ = 4*π*1e-7     # [H m⁻¹] Permeability of free space
Z₀ = sqrt(μ₀/ε₀)  # [Ω] Impedance of free space
c = sqrt(1/ε₀/μ₀) # [m s⁻¹] Speed of light
v = 2/3 * c       # [m s⁻¹] Approximate speed of propagation of lightning strike
                  #        through medium
qₑ = 1.602e-19    # [coul] Charge of electron
mₑ = 9.109e-31    # [kg] Mass of electron
kB = 8.617e-5     # [eV K⁻¹] Boltzmann constant
# Note that these ρ₀ and z₀ are chosen for numerical stability.
# For a rectangular grid discretization, good values are 3Δρ, 3Δz.
#z₀ = 3.0e3        # [m] Heidler normalizing constant
#ρ₀ = 3.0e3        # [m] Heidler normalizing constant
z₀ = 3 * (0.5e3)  # [m] Heidler normalizing constant
ρ₀ = 3 * (0.5e3)    # [m] Heidler normalizing constant
a = 10.0 * 1e3    # [m] Starting altitude of vertical decay for Heidler
# Note: If the axes are swapped, then these accessors must be as well.
ρ = EForm(avg_mat * VForm(map(p -> p[1], sd[:point]))) # [m] Distance from centre of strike
Z = EForm(avg_mat * VForm(map(p -> p[2], sd[:point]))) # [m] Altitude

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
# Heidler Model #
#################
begin # Heidler Model

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
  # See Kotovsky pp. 92 5-2b, and Pasko et al. ...
  Eq5_2b == qₑ * 3.656e25 * ρ_e ./ ρ_gas .* sqrt(200.0 ./ Tn)

  σ == avg₀₁((Eq5_2a_mask .*  Eq5_2a) + (Eq5_2b_mask .*  Eq5_2b))
end
end # Veronis

################################################
# Initialize Densities and Neutral Temperature #
################################################
begin # Format Atmosphere

include("./examples/LightnignChanges copy/formatAtmosphere.jl")
species, densities, Tn, rateCoef_names, rateCoefs = formatAtmosphere("./examples/LightnignChanges copy/chi180_O_dyn.mat", sd)
#mesh(s, color=densities[:O])

end # Format Atmosphere

####################
# Initialize E & B #
####################
begin # Initialize Electromagnetism

# Note: we are using the Linear model of the reduced electric field
E₀ = VForm(zeros(ne(s)))
B₀ = DualForm{0}(zeros(ntriangles(s)))

end # Initialize Electromagnetism

#####################
# Model Composition #
#####################
begin # Model Composition

# TODO: Rename this file to chemistry.jl
include("./examples/LightnignChanges copy/rateCoefficients_dynam.jl")

# This infers ~750 more forms.
@show count(==(:infer), Chemistry[:type])
rate_idxs = findall(x -> String(x)[1] == 'r', Chemistry[:name]);
coef_idxs = findall(x -> String(x)[1] == 'k', Chemistry[:name]);
Chemistry[rate_idxs, :type] = fill(:Form0, length(rate_idxs));
Chemistry[coef_idxs, :type] = fill(:Form0, length(coef_idxs));
@show count(==(:infer), Chemistry[:type])
#infer_types!(Chemistry);

#bespoke_op1_inf_rules = [
#(src_type = :Form0, tgt_type = :Form0, op_names = [:exp])]
#infer_types!(Chemistry, vcat(bespoke_op1_inf_rules, op1_inf_rules_2D),
#  op2_inf_rules_2D)

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
@show count(==(:infer), lightning[:type])

end # Model Composition

#########################
# Simulation Generation #
#########################
begin # Simulation Generation

#open("./generated_lightning_sim_sep21.jl", "w") do file
#  write(file, string(gensim(lightning)))
#end
sim = include("./generated_lightning_sim_sep21.jl")
fₘ = sim(sd, generate, DiagonalHodge())

end

###########
# Solving #
###########

u₀ = construct(PhysicsState,
  map(x -> VectorForm(x.data),
    [Z,
     ρ,
     B₀,
     E₀,
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

# O2, e, O
function check_neg_densities(u)
  for sp in species
    sp_key = Symbol(sp == :e ? :ρ_e : "Chemistry_ρ_" * string(sp))
    ρ = findnode(u, sp_key)
    any(ρ .< 0.0) && println(string(sp_key))
  end
end

du₀ = deepcopy(u₀)
dt = 3e-9
# tₑ = 200e-6
tₑ = 76e-6
fₘ(du₀, u₀, constants_and_parameters, 0)
using Profile
using ProfileSVG
Profile.clear()
@profile fₘ(du₀, u₀, constants_and_parameters, 0)
ProfileSVG.save("prof_lightning_sum.svg")
using BenchmarkTools
@btime fₘ(du₀, u₀, constants_and_parameters, 0)
using ProgressBars
for t in ProgressBar(0:dt:tₑ)
  fₘ(du₀, u₀, constants_and_parameters, t)
  u₀ .+= du₀*dt
  check_neg_densities(u₀)
end

# Chemistry_P_O .= (.+)(Chemistry_r7, ...)
#julia> @btime fₘ(du₀, u₀, constants_and_parameters, 0)
#  762.302 ms (18179595 allocations: 709.49 MiB)
# Chemistry_P_O .= sum([Chemistry_r7, ...])
#julia> @btime fₘ(du₀, u₀, constants_and_parameters, 0)
#  133.040 ms (629 allocations: 93.61 MiB)

# These are typical values for a model run:
#tₑ = 200e-6 # [s]
#tₑ = 800e-6 # [s]

# This is how long it takes light to travel 400 km in a vacuum.
#tₑ = 1.334e-3 # [s] 

tₑ = 200e-6 * .67 # tₑ =  76e-6
tₑ = 200e-6 * .3
#tₑ = 200e-6 * .67 * .57
tₑ = 76e-6
prob = ODEProblem{true}(fₘ, u₀, (0, tₑ), constants_and_parameters)
@time soln = solve(prob, Tsit5(), progress=true, progress_steps=1, save_everystep=false)
# 60e-6
#277.856911 seconds (4.42 G allocations: 212.818 GiB, 6.80% gc time)
@time soln = solve(prob, Tsit5(), progress=true, progress_steps=1)

# 60e-6
#285.874423 seconds (4.42 G allocations: 228.682 GiB, 6.38% gc time, 0.05% compilation time)
#julia> extrema(soln.t[2:end] - soln.t[1:end-1])
#(6.53460176053054e-7, 1.9230405967132077e-5)

#tₑ = 76e-6:
#425.802378 seconds (8.02 G allocations: 414.330 GiB, 5.09%
# gc time, 0.08% compilation time: <1% of which was recompilation)
#tₑ = 76e-6, with pd_wedge for magnitude of E:

#tₑ = 76e-6:
#@time soln = solve(prob, Tsit5(), progress=true, progress_steps=1, save_everystep=false)
#ODE 100%|██████████████████████████████████████████████████████████████████████████████████████| Time: 0:09:27
#573.011931 seconds (8.07 G allocations: 314.896 GiB, 3.03% gc time, 28.26% compilation time: <1% of which was recompilation)

# After using sum instead of (.+)
#tₑ = 76e-6:
#julia> @time soln = solve(prob, Tsit5(), progress=true, progress_steps=1, save_everystep=false)
#ODE 100%|█████████████████████████████████████████████████████████████████| Time: 0:01:25
# 85.400361 seconds (396.79 k allocations: 46.359 GiB, 0.34% gc time)

# And further setting I₀ to 250.
#tₑ = 76e-6:
#
#julia> @time soln = solve(prob, Tsit5(), progress=true, progress_steps=1, save_everystep=false)
#ODE 100%|█████████████████████████████████████████████████████████████████| Time: 0:00:58
# 58.694405 seconds (277.96 k allocations: 32.257 GiB, 0.11% gc time)

tₑ = 3e-9 * 16
prob = ODEProblem{true}(fₘ, u₀, (0, tₑ), constants_and_parameters)
@time soln = solve(prob, ORK256(), dt=3e-9, progress=true, progress_steps=10, save_everystep=false, saveat=[1:10..., 20:10:60..., 61:76...]*1e-6)
@time soln = solve(prob, ORK256(), dt=3e-9, progress=true, progress_steps=1, save_everystep=false)

prob = ODEProblem{true}(fₘ, u₀, (0, tₑ), constants_and_parameters)
@time soln = solve(prob, Tsit5(), dtmax=3e-9, progress=true, progress_steps=1, save_everystep=false)

# Warning: This overwrites a file with this name if there is one.
@save "/orange/fairbanksj/luke.morris/lightning_data/lightning_60em6_dense.jld2" soln
@save "/orange/fairbanksj/luke.morris/lightning_data/lightning_76em6_dense.jld2" soln
@save "/orange/fairbanksj/luke.morris/lightning_data/lightning_76em6_sep21.jld2" soln
@save "/orange/fairbanksj/luke.morris/lightning_data/lightning_76em6_lowIO.jld2" soln

# There is no need to load the solution back in if you haven't killed the REPL.
@load "lightning.jld2" soln 
@load "/orange/fairbanksj/luke.morris/lightning_data/lightning_60em6_dense.jld2" soln
@load "/orange/fairbanksj/luke.morris/lightning_data/lightning_76em6_dense.jld2" soln
tₑ = last(soln.t)

############
# Plotting #
############

# Compute the divergence. i.e. ⋆d⋆ in the Discrete Exterior Calculus
div_mat = inv_hodge_star(0,sd)*dual_derivative(1,sd)*hodge_star(1,sd)
div(form1) = div_mat*form1

# Plot the divergence of the E field. i.e. ⋆d⋆(E)
function plot_div_E(t, colorrange=extrema(findnode(soln(t), :Veronis_E)))
  mesh(s, color=div(findnode(soln(t), :Veronis_E)); colorrange=colorrange)
end

plots_prefix = "/orange/fairbanksj/luke.morris/lightning_data/"
time_string = "76em6"

fig = plot_div_E(tₑ)
save(plots_prefix * "/div_E_" * time_string * ".png", fig)
fig = plot_div_E(tₑ, (-1e-32, 1e-32))
save(plots_prefix * "/div_E_cr_" * time_string * ".png", fig)

# Plot the density of e over the initial density of e.
function plot_px_e(t)
  mesh(s, color=findnode(soln(tₑ), :ρ_e) / findnode(soln(0), :ρ_e))
end

plot_px_e(tₑ)

# Plot the density of O over the initial density of O.
function plot_px_O(t)
  mesh(s, color=findnode(soln(tₑ), :Chemistry_ρ_O) / findnode(soln(0), :Chemistry_ρ_O))
end

plot_px_O(tₑ)

function plot_all_diffs(species_considered = species)
  fig = Figure(fontize=24)
  for ((i,sp),fig_idx) in zip(enumerate(species_considered),
    Iterators.product(Int.(1:sqrt(length(species_considered))+1),Int.(1:sqrt(length(species_considered))+1)))
    sp_key = Symbol(sp == :e ? :ρ_e : "Chemistry_ρ_" * string(sp))
    diff = findnode(soln(tₑ), sp_key) - findnode(soln(0), sp_key)
    ax = Axis(fig[fig_idx...],
              title="Difference in Density of " * String(sp),
              aspect=1,
              width=800, height=800)
    #colsize!(fig.layout, 1, Aspect(1, 4.0))
    mesh!(ax, s, color=diff)
  end
  resize_to_layout!(fig)
  fig
end

fig = plot_all_diffs()
save(plots_prefix * "/all_diffs_" * time_string * ".png", fig)
dynamic_species = [:N2A, :N2B, :N2a, :N2C, :N4S, :N2D, :O2a, :O2b, :NO, :NO2, :e, :Om, :O2m, :O3m, :O4m, :OHm, :CO3m, :CO4m, :NO2m, :NO3m, :O2mNO, :HCO3m, :N2p, :Np, :O2p, :Op, :O4p, :NOp, :Yp, :O]
fig = plot_all_diffs(dynamic_species)
save(plots_prefix * "/all_diffs_w_Ps_Ls_" * time_string * ".png", fig)

function plot_all_diffs_above_60(species_considered = species)
  fig = Figure(fontize=24)
  points_mask = map(x -> x[2], point(s)) .> 60e3
  for ((i,sp),fig_idx) in zip(enumerate(species_considered),
    Iterators.product(Int.(1:sqrt(length(species_considered))+1),Int.(1:sqrt(length(species_considered))+1)))
    sp_key = Symbol(sp == :e ? :ρ_e : "Chemistry_ρ_" * string(sp))
    diff = findnode(soln(tₑ), sp_key) - findnode(soln(0), sp_key)
    ax = Axis(fig[fig_idx...],
              title="Difference in Density of " * String(sp),
              aspect=1,
              width=800, height=800)
    scatter!(ax, point(s)[points_mask], color=diff[points_mask])
  end
  resize_to_layout!(fig)
  fig
end

fig = plot_all_diffs_above_60()
save(plots_prefix * "/all_diffs_above_60_" * time_string * ".png", fig)
fig = plot_all_diffs_above_60([:N2A, :N2B, :N2a, :N2C, :N4S, :N2D, :O2a, :O2b, :NO, :NO2, :e, :Om, :O2m, :O3m, :O4m, :OHm, :CO3m, :CO4m, :NO2m, :NO3m, :O2mNO, :HCO3m, :N2p, :Np, :O2p, :Op, :O4p, :NOp, :Yp, :O])
save(plots_prefix * "/all_diffs_w_Ps_Ls_above_60_" * time_string * ".png", fig)

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

# Plot the E field itself.
function plot_E(t)
  fig = Figure(resolution=(1e3,1e3))
  ax = Axis(fig[1,1])
  max_x = find_max_x_nonzero_form1(findnode(soln(t), :Veronis_E))
  points_mask = map(x -> x[1], point(s)) .< max_x
  CairoMakie.arrows!(ax,
                     point(s)[points_mask],
                     (sharp_mat*findnode(soln(t), :Veronis_E))[points_mask],
                     lengthscale=1e-6)
  fig
end

fig = plot_E(tₑ)
save(plots_prefix * "/E_arrows" * time_string * ".png", fig)

# Plot the magnitude of the E field.
function plot_mag_E(t)
  mesh(s, color=norm.(sharp_mat*findnode(soln(t), :Veronis_E)))
end

fig = plot_mag_E(tₑ)
save(plots_prefix * "/E_mag" * time_string * ".png", fig)

# Plot the log of the magnitude of the E field.
function plot_log_mag_E(t)
  f = Figure()
  ax = Axis(f[1,1],
            title = "ln||E||",
            yticks = 0:(MAX_Z/10):MAX_Z)
  mesh!(ax, s, color=log.(norm.(sharp_mat*findnode(soln(t), :Veronis_E))))
  f
end

fig = plot_log_mag_E(tₑ)
save(plots_prefix * "/E_log_mag" * time_string * ".png", fig)

# Plot the log of the magnitude of the E field, bounded.
function plot_log_mag_E_bounded(t)
  f = Figure()
  points_mask = map(x -> x[1], point(s)) .< 150e3
  ax = Axis(f[1,1],
            title = "ln||E||",
            yticks = 0:(MAX_Z/10):MAX_Z)
  sctr = scatter!(ax, point(s)[points_mask], color=(log.(norm.(sharp_mat*findnode(soln(t), :Veronis_E))))[points_mask])
  Colorbar(fig[1,2], sctr)
  f
end

fig = plot_log_mag_E_bounded(tₑ)
save(plots_prefix * "/E_log_mag_bounded_" * time_string * ".png", fig)

# Plot the log of the magnitude of J.
function plot_log_mag_J(t)
  f = Figure()
  ax = Axis(f[1,1],
            title = "ln||J||",
            yticks = 0:(MAX_Z/10):MAX_Z)
  mesh!(ax, s, color=log.(norm.(sharp_mat*J)))
  f
end

fig = plot_log_mag_J(tₑ)
save(plots_prefix * "/J_log_mag" * time_string * ".png", fig)

# Plot the log of the absolute value of the B field.
function plot_log_abs_B(t)
  f = Figure()
  ax = Axis(f[1,1],
            title = "ln|B|",
            yticks = 0:(MAX_Z/10):MAX_Z)
  scatter!(ax, dual_point(sd, triangle_center(sd, triangles(sd))),
    color=(log ∘ abs).(findnode(soln(t), :Veronis_B)))
  f
end

fig = plot_log_abs_B(tₑ)
save(plots_prefix * "/B_log_abs" * time_string * ".png", fig)

# Plot the B field itself.
function plot_B(t)
  scatter(dual_point(sd, triangle_center(sd, triangles(sd))),
    color=findnode(soln(t), :Veronis_B))
end

plot_B(tₑ)

# Create a gif of magnitude of E.
# This requires that you did not set save_everystep to false when you solved.
function make_gif_E()
  frames = 100
  fig, ax, ob = mesh(s, color=norm.(sharp_mat * findnode(soln(0), :Veronis_E))), colormap=:jet, colorrange=extrema(norm.(sharp_mat * findnode(soln(tₑ), :Veronis_E)))
  Colorbar(fig[1,2], ob)
  record(fig, "lightning.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    ob.color = norm.(♯_m * findnode(soln(t), :Veronis_E))
  end
end

make_gif_E()
  fig, ax, ob = mesh(s, color=norm.(sharp_mat * findnode(soln(0), :Veronis_E))), colormap=:jet, colorrange=extrema(norm.(sharp_mat * findnode(soln(tₑ), :Veronis_E)))
  Colorbar(fig[1,2], ob)
  record(fig, "lightning.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    ob.color = norm.(♯_m * findnode(soln(t), :Veronis_E))
  end
end

make_gif_E()

# Plot 10 * log of the density of the species at time t.
function plot_10_log_ρ_above_60(sp, t)
  sp_key = Symbol(sp == :e ? :ρ_e : "Chemistry_ρ_" * string(sp))
  points_mask = map(x -> x[2], point(s)) .> 60e3
  fig, ax, ob = scatter(point(s)[points_mask], color=10*log10.(findnode(soln(0), sp_key)[points_mask]), colormap=:jet, colorrange=extrema(10*log10.(findnode(soln(0), sp_key)[points_mask])))
  Colorbar(fig[1,2], ob)
  fig
end

fig = plot_10_log_ρ_above_60(:e, tₑ)
save(plots_prefix * "10log10_density_e_above_60_" * time_string * ".png", fig)

# Create a gif of ρ of the given species.
# This requires that you did not set save_everystep to false when you solved.
function make_gif_10_log_ρ(sp, microseconds)
  sp_key = Symbol(sp == :e ? :ρ_e : "Chemistry_ρ_" * string(sp))
  fig, ax, ob = mesh(s, color=10*log10.(findnode(soln(0), sp_key)), colormap=:jet, colorrange=extrema(findnode(soln(0), sp_key)))
  Colorbar(fig[1,2], ob)
  record(fig, plots_prefix * "/10log10_density_" * String(sp) * "_" * time_string * ".gif", microseconds; framerate = 2) do t
    ob.color = 10*log10.(findnode(soln(t*1e-6), sp_key))
    ax.title = "10*log₁₀ of Density of " * String(sp) * " at t = " * string(t)
  end
end

make_gif_10_log_ρ(:e, [1:10..., 20:10:60...])
make_gif_10_log_ρ(:e, [1:10..., 20:10:60..., 61:76...])

# Create a gif of ρ of the given species above 60.
# This requires that you did not set save_everystep to false when you solved.
function make_gif_10_log_ρ_above_60(sp, microseconds)
  sp_key = Symbol(sp == :e ? :ρ_e : "Chemistry_ρ_" * string(sp))
  points_mask = map(x -> x[2], point(s)) .> 60e3
  fig, ax, ob = scatter(point(s)[points_mask], color=10*log10.(findnode(soln(0), sp_key)[points_mask]), colormap=:jet, colorrange=extrema(10*log10.(findnode(soln(0), sp_key)[points_mask])))
  Colorbar(fig[1,2], ob)
  record(fig, plots_prefix * "/10log10_density_" * String(sp) * "_" * time_string * ".gif", microseconds; framerate = 2) do t
    ob.color = 10*log10.(findnode(soln(t*1e-6), sp_key)[points_mask])
    ax.title = "10*log₁₀ of Density of " * String(sp) * " at t = " * string(t) * "usec"
  end
end

make_gif_10_log_ρ_above_60(:e, [1:10..., 20:10:60...])
make_gif_10_log_ρ_above_60(:e, [1:10..., 20:10:60..., 61:76...])

# Heidler waveform plots:
plot(0.0:1.0e-8:100.0e-6, HeidlerForη)
log_Heidler = map(0.0:1.0e-5:1.5e-3) do x
  log(HeidlerForη(x)*1e6)
end
log_Heidler[1] = log_Heidler[2] # Set this manually to something not -Inf
plot(0.0:1.0e-5:1.5e-3, log_Heidler)
plot(0.0:1.0e-8:100.0e-6, ddtHeidlerForη)

# Magnitude of E at a point comparisons.
Ef = findnode(soln(tₑ), :Veronis_E)
invstar0_mat = inv_hodge_star(0, sd)

Ef = E₀
Efs = norm.(sharp_mat * Ef)
fig = mesh(s, color=Efs)
save(plots_prefix * "sharp_mag_E" * time_string * ".png", fig)
fig = mesh(s, color=log10.(Efs))
save(plots_prefix * "log_sharp_mag_E" * time_string * ".png", fig)

Efw = invstar0_mat * pd_wedge(Val{(1,1)}, sd, Ef, star1_mat * Ef)
fig = mesh(s, color=Efw)
save(plots_prefix * "wedge_mag_E" * time_string * ".png", fig)
fig = mesh(s, color=log10.(Efw))
save(plots_prefix * "log_wedge_mag_E" * time_string * ".png", fig)

#julia> onedx_plus_onedy = eval_constant_form(sd, @SVector [1.0, 1.0, 0.0]);

#julia> extrema(invstar0_mat * pd_wedge(Val{(1,1)}, sd, onedx_plus_onedy, star1_mat * onedx_plus_onedy; wedge_t=wedge_m))
#(1.8325664743647492, 2.167433525635252)

#julia> histogram(invstar0_mat * pd_wedge(Val{(1,1)}, sd, onedx_plus_onedy, star1_mat * onedx_plus_o
#nedy; wedge_t=wedge_m))
#                ┌                                        ┐ 
#   [1.82, 1.84) ┤▏ 1                                       
#   [1.84, 1.86) ┤  0                                       
#   [1.86, 1.88) ┤  0                                       
#   [1.88, 1.9 ) ┤  0                                       
#   [1.9 , 1.92) ┤▏ 1                                       
#   [1.92, 1.94) ┤  0                                       
#   [1.94, 1.96) ┤  0                                       
#   [1.96, 1.98) ┤  0                                       
#   [1.98, 2.0 ) ┤████████████████████████████████  22 343  
#   [2.0 , 2.02) ┤██████████████████████████▌ 18 504        
#   [2.02, 2.04) ┤  0                                       
#   [2.04, 2.06) ┤  0                                       
#   [2.06, 2.08) ┤  0                                       
#   [2.08, 2.1 ) ┤▏ 1                                       
#   [2.1 , 2.12) ┤  0                                       
#   [2.12, 2.14) ┤  0                                       
#   [2.14, 2.16) ┤  0                                       
#   [2.16, 2.18) ┤▏ 1                                       
                └                                        ┘ 
#
#julia> extrema(norm.(sharp_mat * eval_constant_form(sd, @SVector [1.0,1.0,0.0])))
#(1.414213562373094, 1.4142135623730963)
#
#julia> histogram(norm.(sharp_mat * eval_constant_form(sd, @SVector [1.0,1.0,0.0])))
#                                            ┌                                        ┐ 
#                        [1.41421 , 1.41421) ┤▎ 76                                      
#                         [1.41421, 1.41421) ┤▌ 161                                     
#                         [1.41421, 1.41421) ┤██████▍ 2 235                             
#   [0.0               , 0.0               ) ┤█████████▋ 3 362                          
#                        [1.41421, 1.41421 ) ┤████████████████████████▎ 8 448           
#                        [1.41421 , 1.41421) ┤  0                                       
#                         [1.41421, 1.41421) ┤████████████████████████████████  11 175  
#                         [1.41421, 1.41421) ┤███████████████████▌ 6 798                
#                         [1.41421, 1.41421) ┤███████████████████████▏ 8 053            
#                        [1.41421, 1.41421 ) ┤█▍ 452                                    
#                        [1.41421 , 1.41421) ┤▎ 71                                      
#                         [1.41421, 1.41421) ┤▏ 20                                      
#                                            └                                        ┘ 
#                                                             Frequency                 
#julia> extrema(norm.(sharp_mat * onedx_plus_onedy)) .- √2
#(-1.1102230246251565e-15, 1.1102230246251565e-15)
