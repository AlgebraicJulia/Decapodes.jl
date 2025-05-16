# Author: Luke Morris
# This is a discretization of the incompressible Navier Stokes equations using the Discrete Exterior Calculus.
#
# These formulations are based on those given by Mohamed, Hirani, Samtaney, (in turn from Marsden, Ratiu, Abraham.)
#
# However, different choices in discretization are chosen for purposes of brevity, to demonstrate novel discretizations of certain operators, and to demonstrate the automated Decapodes workflow.

################
# Dependencies #
################
begin # Dependencies
@info "Loading Dependencies"

# AlgebraicJulia:
using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations

# Meshing:
using CoordRefSystems
using GeometryBasics: Point3
Point3D = Point3{Float64};

# Visualization:
using CairoMakie

# Simulation:
using ComponentArrays
using LinearAlgebra
using LinearAlgebra: factorize
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using StaticArrays

# Saving:
using JLD2
using Parameters
using TOML

# Logging:
using LoggingExtras
using Logging: global_logger
using Statistics: mean
using TerminalLoggers: TerminalLogger
end # Dependencies

#########
# Flags #
#########
begin # Flags
@info "Setting Flags"
config = TOML.parsefile("config.toml")
@unpack USE_EQ11,USE_EQ11_INVISCID, USE_EQ11_POISSON, USE_EQ17, USE_EQ17_INVISCID = config
@assert xor(USE_EQ11, USE_EQ11_INVISCID, USE_EQ11_POISSON, USE_EQ17, USE_EQ17_INVISCID)

@unpack TAYLOR_SCENARIO, SIX_VORTEX_SCENARIO, CUSTOM_SCENARIO = config
@assert xor(TAYLOR_SCENARIO, SIX_VORTEX_SCENARIO, CUSTOM_SCENARIO)

@unpack CUSTOM_LAT, CUSTOM_N_VORTS, CUSTOM_tau, CUSTOM_G, CUSTOM_a,CUSTOM_POINT_VORTEX, CUSTOM_TAYLOR_VORTEX, CUSTOM_TIME, CUSTOM_CLRRNG = config

@unpack ICO6, ICO7, ICO8, UV = config
@assert xor(ICO6, ICO7, ICO8, UV)

@unpack BARYCENTRIC, CIRCUMCENTRIC = config
@assert xor(BARYCENTRIC, CIRCUMCENTRIC)

@unpack DATA_DR, SAVE_SOLN, PLOT_ICS, PLOT_FCS, SCP_ICS, SCP_FCS, SERVER_PATH = config

@unpack LOG_TERMINAL, LOG_FILE = config
@assert xor(LOG_TERMINAL, LOG_FILE)

if LOG_TERMINAL
  global_logger(TerminalLogger())
elseif LOG_FILE
  logger = FileLogger(joinpath(DATA_DR, "vort.log"); append=true)
  global_logger(logger)
end
end # Flags

##########
# Models #
##########
begin # Models
@info "Defining Models"
eq11_vorticity = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1
  μ::Constant

  𝐮 == d₁⁻¹(d𝐮)

  ∂ₜ(d𝐮) == μ * ∘(⋆, d, ⋆, d)(d𝐮) + (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
open("eq11_vorticity.jl", "w") do f
  write(f, string(gensim(eq11_vorticity)))
end

eq11_inviscid_vorticity = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1

  𝐮 == d₁⁻¹(d𝐮)

  ∂ₜ(d𝐮) ==  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
open("eq11_inviscid_vorticity.jl", "w") do f
  write(f, string(gensim(eq11_inviscid_vorticity)))
end

eq11_inviscid_poisson = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1
  ψ::Form0

  ψ == Δ⁻¹(⋆(d𝐮))
  𝐮 == ⋆(d(ψ))

  ∂ₜ(d𝐮) ==  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end
open("eq11_inviscid_poisson.jl", "w") do f
  write(f, string(gensim(eq11_inviscid_poisson)))
end

eq17_stream = @decapode begin
  ψ::Form0
  u::DualForm1
  v::Form1
  μ::Constant

  u == ⋆(d(ψ))
  v == ⋆(u)

  ∂ₜ(ψ) == dsdinv(
                  μ * ∘(d, ⋆, d, ⋆, d, ⋆, d)(ψ) -
                  ∘(⋆₁, d̃₁)(v ∧ ∘(d,⋆,d,⋆)(ψ)))
end
open("eq17_stream.jl", "w") do f
  write(f, string(gensim(eq17_stream)))
end

eq17_inviscid_stream = @decapode begin
  ψ::Form0
  u::DualForm1
  v::Form1

  u == ⋆(d(ψ))
  v == ⋆(u)

  ∂ₜ(ψ) == -1 * dsdinv(∘(⋆₁, d̃₁)(v ∧ ∘(d,⋆,d,⋆)(ψ)))
end
open("eq17_inviscid_stream.jl", "w") do f
  write(f, string(gensim(eq17_inviscid_stream)))
end

sim = if USE_EQ11
  include("eq11_vorticity.jl")
elseif USE_EQ11_INVISCID
  include("eq11_inviscid_vorticity.jl")
elseif USE_EQ11_POISSON
  include("eq11_inviscid_poisson.jl")
elseif USE_EQ17
  include("eq17_stream.jl")
elseif USE_EQ17_INVISCID
  include("eq17_inviscid_stream.jl")
end
end # Models

######################
# Mesh and Operators #
######################
begin # Mesh & Ops
@info "Allocating Mesh and Operators"
const RADIUS = 1.0
s = if ICO6
  loadmesh(Icosphere(6, RADIUS));
elseif ICO7
  loadmesh(Icosphere(7, RADIUS));
elseif ICO8
  loadmesh(Icosphere(8, RADIUS));
elseif UV
  s, _, _ = makeSphere(0, 180, 2.5, 0, 360, 2.5, RADIUS);
  s;
end;
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
if BARYCENTRIC
  subdivide_duals!(sd, Barycenter());
elseif CIRCUMCENTRIC
  subdivide_duals!(sd, Circumcenter());
end
@info "    Mesh of $(ntriangles(sd)) triangles allocated"

d0 = dec_differential(0,sd);
d1 = dec_differential(1,sd)
dd0 = dec_dual_derivative(0,sd);
dd1 = dec_dual_derivative(1,sd);
fd0 = factorize(float.(d0));
fdd1 = factorize(float.(dd1));
δ1 = δ(1,sd)
s0 = dec_hodge_star(0,sd,GeometricHodge())
s1 = dec_hodge_star(1,sd,GeometricHodge())
s2 = dec_hodge_star(2, sd)
s0inv = dec_inv_hodge_star(0,sd,GeometricHodge())
Δ0 = Δ(0,sd);
fΔ0 = factorize(Δ0);
♭♯_m = ♭♯_mat(sd);
dsd = factorize(dd1 * s1 * d0);
# As defined in the MHS paper:
dᵦ = 0.5 * abs.(dd1) * spdiagm(dd0 * ones(ntriangles(sd)));
@info "    Differential operators allocated"

function generate(s, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :d₁⁻¹ => x -> fdd1 \ x
    :Δ⁻¹ => x -> begin
      y = fΔ0 \ x
      y .- minimum(y)
    end
    :dsdinv => x -> dsd \ x
    :dinv => x -> fd0 \ x
    :♭♯ => x -> ♭♯_m * x
    :dᵦ => x -> dᵦ * x
    _ => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end;

end

#######################
# Generate Simulation #
#######################
begin # Generate Simulation
@info "Generating Simulation"
fₘ = sim(sd, generate);
end # Generate Simulation

######################
# Initial Conditions #
######################
begin # ICs
@info "Setting Initial Conditions"

"""    function great_circle_dist(pnt,G,a,cntr)

Compute the length of the shortest path along a sphere, given Cartesian coordinates.
"""
function great_circle_dist(pnt1::Point3D, pnt2::Point3D)
  RADIUS * acos(dot(pnt1,pnt2))
end

abstract type AbstractVortexParams end

struct TaylorVortexParams <: AbstractVortexParams
  G::Real
  a::Real
end

struct PointVortexParams <: AbstractVortexParams
  τ::Real
  a::Real
end

"""    function taylor_vortex(pnt::Point3D, cntr::Point3D, p::TaylorVortexParams)

Compute the value of a Taylor vortex at the given point.
"""
function taylor_vortex(pnt::Point3D, cntr::Point3D, p::TaylorVortexParams)
  gcd = great_circle_dist(pnt,cntr)
  (p.G/p.a) * (2 - (gcd/p.a)^2) * exp(0.5 * (1 - (gcd/p.a)^2))
end

"""    function point_vortex(pnt::Point3D, cntr::Point3D, p::PointVortexParams)

Compute the value of a smoothed point vortex at the given point.
"""
function point_vortex(pnt::Point3D, cntr::Point3D, p::PointVortexParams)
  gcd = great_circle_dist(pnt,cntr)
  p.τ / (cosh(3gcd/p.a)^2)
end

taylor_vortex(sd::HasDeltaSet, cntr::Point3D, p::TaylorVortexParams) =
  map(x -> taylor_vortex(x, cntr, p), point(sd))
point_vortex(sd::HasDeltaSet, cntr::Point3D, p::PointVortexParams) =
  map(x -> point_vortex(x, cntr, p), point(sd))

"""    function ring_centers(lat, n)

Find n equispaced points at the given latitude.
"""
function ring_centers(lat, n)
  ϕs = range(0.0, 2π; length=n+1)[1:n]
  map(ϕs) do ϕ
    v_sph = Spherical(RADIUS, lat, ϕ)
    v_crt = convert(Cartesian, v_sph)
    Point3D(v_crt.x.val, v_crt.y.val, v_crt.z.val)
  end
end

"""    function vort_ring(lat, n_vorts, p::T, formula) where {T<:AbstractVortexParams}

Compute vorticity as primal 0-forms for a ring of vortices.

Specify the latitude, number of vortices, and a formula for computing vortex strength centered at a point.
"""
function vort_ring(lat, n_vorts, p::T, formula) where {T<:AbstractVortexParams}
  sum(map(x -> formula(sd, x, p), ring_centers(lat, n_vorts)))
end

"""    function vort_ring(lat, n_vorts, p::PointVortexParams, formula)

Compute vorticity as primal 0-forms for a ring of vortices.

Specify the latitude, number of vortices, and a formula for computing vortex strength centered at a point.

Additionally, place a counter-balance vortex at the South Pole such that the integral of vorticity is 0.
"""
function vort_ring(lat, n_vorts, p::PointVortexParams, formula)
  Xs = sum(map(x -> formula(sd, x, p), ring_centers(lat, n_vorts)))
  Xsp = point_vortex(sd, Point3D(0.0, 0.0, -1.0), PointVortexParams(-1*n_vorts*p.τ, p.a))
  Xs + Xsp
end

X = if TAYLOR_SCENARIO
  # "The centers of the two vortices are separated by a distance
  # of 0.4."
  # Recall that θ is the arc-distance to the north pole.
  vort_ring(0.2, 2, TaylorVortexParams(0.5, 0.1), taylor_vortex)
elseif SIX_VORTEX_SCENARIO
  # Six equidistant points at latitude θ=0.4.
  # "... an additional vortex, with strength τ=-18 and a radius a=0.15, is
  # placed at the south pole (θ=π)."
  vort_ring(0.4, 6, PointVortexParams(3.0, 0.15), point_vortex)
elseif CUSTOM_SCENARIO
  if CUSTOM_POINT_VORTEX
    vort_ring(CUSTOM_LAT, CUSTOM_N_VORTS, PointVortexParamsVortexParams(CUSTOM_tau, CUSTOM_a), point_vortex)
  elseif CUSTOM_TAYLOR_VORTEX
    vort_ring(CUSTOM_LAT, CUSTOM_N_VORTS, TaylorVortexParamsVortexParams(CUSTOM_G, CUSTOM_a), taylor_vortex)
  end
end

"""    function solve_poisson(vort::VForm)

Compute the stream function by solving the Poisson equation.
"""
function solve_poisson(vort::VForm)
  ψ = fΔ0 \ vort.data
  ψ = ψ .- minimum(ψ)
end
solve_poisson(vort::CombinatorialSpaces.DualForm{2}) =
  solve_poisson(VForm(s0inv * vort.data))

ψ = solve_poisson(VForm(X))

# Compute velocity as curl (⋆d) of the stream function.
curl_stream(ψ) = s1 * d0 * ψ
div(u) = s2 * d1 * (s1 \ u)
RMS(x) = √(mean(x' * x))

integral_of_curl(curl::CombinatorialSpaces.DualForm{2}) = sum(curl.data)
# Recall that s0 effectively multiplies each entry by a solid angle.
# i.e. (sum ∘ ⋆₀) computes a Riemann sum.
integral_of_curl(curl::VForm) = integral_of_curl(CombinatorialSpaces.DualForm{2}(s0*curl.data))

u₀ = if USE_EQ11
  # d𝐮::DualForm2
  ComponentArray(d𝐮 = s0*X)
elseif USE_EQ11_INVISCID
  # d𝐮::DualForm2
  ComponentArray(d𝐮 = s0*X)
elseif USE_EQ11_POISSON
  # d𝐮::DualForm2
  ComponentArray(d𝐮 = s0*X)
elseif USE_EQ17
  # ψ::Form0
  ComponentArray(ψ = ψ)
elseif USE_EQ17_INVISCID
  # ψ::Form0
  ComponentArray(ψ = ψ)
end

constants_and_parameters = (
  μ = 0.001,)

@info "RMS of divergence of initial velocity: $(∘(RMS, div, curl_stream)(ψ))"
@info "Integral of initial curl: $(integral_of_curl(VForm(X)))"
end # ICs

###############################
# Plotting Initial Conditions #
###############################
begin # Plot ICs
if PLOT_ICS
  @info "Plotting Initial Conditions"
  function save_ics(ics, ics_name, ics_png_name)
    fig = Figure()
    Label(fig[1, 1, Top()], "$ics_name at t=0", padding = (0, 0, 5, 0))
    ax = LScene(fig[1,1], scenekw=(lights=[],))
    msh = CairoMakie.mesh!(ax, s,
      color=ics,
      colormap=Reverse(:redsblues))
    Colorbar(fig[1,2], msh)
    save(joinpath(DATA_DR, ics_png_name), fig)
  end
  vort_png_name = "vort_ics.png"
  save_ics(X, "Vorticity, X,", vort_png_name)
  if SCP_ICS
    @info "SCPing Initial Conditions Plot"
    run(`scp $DATA_DR/$vort_png_name $SERVER_PATH`)
    @info "    Plot available at $SERVER_PATH/$vort_png_name"
  end
end
end # Plot ICs

############
# Simulate #
############
begin # Simulate
@info("Solving")
tₑ = if TAYLOR_SCENARIO
  10.0
elseif SIX_VORTEX_SCENARIO
  12.0
elseif CUSTOM_SCENARIO
  10.0
end

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob,
  Tsit5(),
  dtmax = 0.01,
  saveat=tₑ/10.0,
  dense=false,
  progress=true, progress_steps=1);
@show soln.retcode
@info("Done")

rms_init = ∘(RMS, div, curl_stream, solve_poisson)(CombinatorialSpaces.DualForm{2}(soln(0).d𝐮))
rms_final = ∘(RMS, div, curl_stream, solve_poisson)(CombinatorialSpaces.DualForm{2}(soln(tₑ).d𝐮))
rms_reldiff = (rms_final - rms_init) / rms_init
@info "RMS of divergence of initial velocity: $rms_init"
@info "RMS of divergence of final velocity: $rms_final"
@info "Relative difference of RMS of divergence: $rms_reldiff"

curl_init = integral_of_curl(CombinatorialSpaces.DualForm{2}(soln(0).d𝐮))
curl_final = integral_of_curl(CombinatorialSpaces.DualForm{2}(soln(tₑ).d𝐮))
curl_reldiff = (curl_final - curl_init) / curl_init
@info "Integral of initial curl: $curl_init"
@info "Integral of final curl: $curl_final"
@info "Relative difference of curl: $curl_reldiff"

if SAVE_SOLN
  @info "Saving Simulation"
  @save joinpath(DATA_DR, "ns.jld2") soln
end
end # Simulate

#########################
# Plot Final Conditions #
#########################
begin # Plot FCs
if PLOT_FCS
  @info "Plotting Final Conditions"

  """    function save_vort_gif(file_name)

  Given a solution wih vorticity as a dual 2-form, make a GIF with vorticity as a primal 0-form.
  """
  function save_vort_gif(file_name, soln)
    time = Observable(0.0)
    fig = Figure()
    Label(fig[1, 1, Top()], @lift("Vorticity at $($time)"), padding = (0, 0, 5, 0))
    ax = LScene(fig[1,1], scenekw=(lights=[],))
    clrrng = if TAYLOR_SCENARIO
      (-5.0,15.0)
    elseif SIX_VORTEX_SCENARIO
      (0.0,3.0)
    elseif CUSTOM_SCENARIO
      CUSTOM_CLRRNG
    end
    msh = CairoMakie.mesh!(ax, s,
      # Observe that s0inv converts d𝐮 from a dual 2-form to a primal 0-form.
      color=@lift(s0inv*soln($time).d𝐮),
      colorrange=clrrng,
      colormap=Reverse(:redsblues))

    Colorbar(fig[1,2], msh)
    record(fig, file_name, soln.t; framerate = 10) do t
      time[] = t
    end
  end

  """    function plot_vort_fc()

  Given a solution wih vorticity as a dual 2-form, plot the final conditions of vorticity as a primal 0-form.
  """
  function plot_vort_fc()
    fig = Figure(size=(2000,2000), fontsize=64)
    Label(fig[1, 1, Top()], "Vorticity at $(last(soln.t))", padding = (0, 0, 5, 0))
    ax = LScene(fig[1,1], scenekw=(lights=[],), show_axis=false)
    update_cam!(ax.scene, Vec3f(0,0,0.8), Vec3f(0,0,0), Vec3f(0, 1, 1))
    clrrng = if TAYLOR_SCENARIO
      (-5.0,15.0)
    elseif SIX_VORTEX_SCENARIO
      (0.0,3.0)
    elseif CUSTOM_SCENARIO
      CUSTOM_CLRRNG
    end
    msh = CairoMakie.mesh!(ax, s,
      # Observe that s0inv converts d𝐮 from a dual 2-form to a primal 0-form.
      color=s0inv*soln(last(soln.t)).d𝐮,
      colorrange=clrrng,
      colormap=Reverse(:redsblues))

    Colorbar(fig[1,2], msh, size=32)
    fig
  end

  """    function azimuth_form(lat, tol, form_ic, form_fc, file_name)

  Plot the vorticity of points lying within the tolerance of the given latitude.

  Give form_ic and form_fc as primal 0-forms.
  """
  function azimuth_form(lat, tol, form_ic::VForm, form_fc::VForm)
    sph_pnts = map(point(sd)) do pnt
      convert(Spherical, Cartesian(pnt...))
    end
    in_lat_range = findall(sph_pnts) do pnt
      abs(pnt.θ - lat) ≤ tol
    end
    by_phi = sortperm(sph_pnts[in_lat_range], by=(x -> x.ϕ))

    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1],
      title="Vorticity along θ=0.4",
      xlabel="Azimuth Angle, φ")
    lns_ic = CairoMakie.lines!(ax,
                               range(0.0,2π; length=length(in_lat_range)),
                               form_ic.data[in_lat_range][by_phi])
    lns_fc = CairoMakie.lines!(ax,
                               range(0.0,2π; length=length(in_lat_range)),
                               form_fc.data[in_lat_range][by_phi])
    Legend(fig[1,2], [lns_ic, lns_fc], ["T=0.0", "T=$(tₑ)"])
    fig
  end

  """    function plot_diff_from_orig()

  Plot the relative difference between the current vorticity and initial vorticity through the length of the simulation.
  """
  function plot_diff_from_orig()
    diff_orig = map(range(0.0, tₑ; length=64)) do t
      norm(soln(t).d𝐮 .- soln(0).d𝐮) / norm(soln(0).d𝐮)
    end
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1],
      title="Relative Solution Change",
      xlabel="Time")
    CairoMakie.lines!(ax, diff_orig)
    fig
  end

  gif_name = "vort.gif"
  save_vort_gif(joinpath(DATA_DR, gif_name), soln)

  fc_plot_name = "vort_final.png"
  save(joinpath(DATA_DR, fc_plot_name),
    plot_vort_fc())

  azimuth_name = "azimuth.png"
  save(joinpath(DATA_DR, azimuth_name),
    azimuth_form(0.4, 0.01, VForm(s0inv*soln(0).d𝐮), VForm(s0inv*soln(tₑ).d𝐮)))

  diff_name = "relsolchng.png"
  save(joinpath(DATA_DR, diff_name),
    plot_diff_from_orig())

  if SCP_FCS
    @info "SCPing Final Conditions"
    foreach([gif_name, fc_plot_name, azimuth_name, diff_name]) do file_name
      run(`scp $DATA_DR/$file_name $SERVER_PATH`)
      @info "    Visualization available at $SERVER_PATH/$file_name"
    end
  end
end
end # Plot FCs
