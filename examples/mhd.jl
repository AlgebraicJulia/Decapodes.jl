# "Magnetohydrodynamics Simulation via Discrete Exterior Calculus", Gillespie, M.
# https://markjgillespie.com/Research/MHD/MHD_Simulation_with_DEC.pdf
# modeled by Matt Cuffaro, Luke Morris
@info "Loading Dependencies"

# AlgebraicJulia
using Catlab
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes

# Meshing
using CoordRefSystems
using GeometryBasics: Point3
Point3D = Point3{Float64};

# Visualization
using CairoMakie

# Simulation
using ComponentArrays
using LinearAlgebra
using LinearAlgebra: factorize
using OrdinaryDiffEq
using SparseArrays
using StaticArrays

# Saving

# other dependencies
using MLStyle
using Statistics: mean

@info "Defining models"
_mhd = @decapode begin
    ψ::Form0
    η::DualForm1
    (dη,β)::DualForm2
    # δ = ⋆d⋆
    ∂ₜ(dη) == -1*(∘(⋆₁, dual_d₁)((⋆(dη) ∧₀₁ ♭♯(η)) + (⋆(β) ∧₀₁ ♭♯(∘(⋆, d, ⋆)(β)))))
    ∂ₜ(β) == -1*(∘(⋆₁, dual_d₁)(⋆(β) ∧₀₁ ♭♯(η)))
    # solve for stream function
    ψ == ∘(⋆, Δ⁻¹)(dη)
    η == ⋆(d(ψ))
end;

mhd = @decapode begin
    ψ::Form0
    (η,β)::DualForm1
    dη::DualForm2
    # δ = ⋆d⋆ # TODO Luke make primal-dual 1, 0
    ∂ₜ(dη) == -1*(∘(⋆₁, dual_d₁)((⋆(dη) ∧₀₁ ♭♯(η)) + (♭♯(⋆(β)) ∧ᵈᵈ₁₀ ∘(⋆, d, ⋆)(β))))
    ∂ₜ(β) == -1*(∘(⋆₂, dual_d₀)(⋆(β) ∧ᵖᵈ₁₁ η))
    # solve for stream function
    ψ == ∘(⋆, Δ⁻¹)(dη)
    η == ⋆(d(ψ))
end;

@info "Allocating Mesh and Operators"
const RADIUS = 1.0;
sphere = :ICO7;
s = @match sphere begin
    :ICO6 => loadmesh(Icosphere(6, RADIUS));
    :ICO7 => loadmesh(Icosphere(7, RADIUS));
    :ICO8 => loadmesh(Icosphere(8, RADIUS));
    :flat => triangulated_grid(10, 10, 0.2, 0.2, Point3D)
    :UV => begin
        s, _, _ = makeSphere(0, 180, 2.5, 0, 360, 2.5, RADIUS);
        s;
    end
end;
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

d0 = dec_differential(0,sd);
d1 = dec_differential(1,sd);
dd0 = dec_dual_derivative(0,sd);
dd1 = dec_dual_derivative(1,sd);
fd0 = factorize(float.(d0));
fdd1 = factorize(float.(dd1));
δ1 = δ(1,sd);
s0 = dec_hodge_star(0,sd,GeometricHodge());
s1 = dec_hodge_star(1,sd,GeometricHodge());
s2 = dec_hodge_star(2, sd);
s0inv = dec_inv_hodge_star(0,sd,GeometricHodge());
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
    _ => default_dec_matrix_generate(s, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end;

sim = gensim(mhd);
open("mhd_sim.jl", "w") do f
    write(f, string(gensim(mhd)))
end
sim = include("mhd_sim.jl")
f = sim(sd, generate);

constants_and_parameters = (
  μ = 0.001,)

###################
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

X =  # Six equidistant points at latitude θ=0.4.
  # "... an additional vortex, with strength τ=-18 and a radius a=0.15, is
  # placed at the south pole (θ=π)."
  vort_ring(0.4, 6, PointVortexParams(3.0, 0.15), point_vortex)

"""    function solve_poisson(vort::VForm)
Compute the stream function by solving the Poisson equation.
"""
function solve_poisson(vort::VForm)
  ψ = fΔ0 \ vort.data
  ψ = ψ .- minimum(ψ)
end
solve_poisson(vort::DualForm{2}) =
  solve_poisson(VForm(s0inv * vort.data))

ψ = solve_poisson(VForm(X))

# Compute velocity as curl (⋆d) of the stream function.
curl_stream(ψ) = s1 * d0 * ψ
div(u) = s2 * d1 * (s1 \ u)
RMS(x) = √(mean(x' * x))

integral_of_curl(curl::DualForm{2}) = sum(curl.data)
# Recall that s0 effectively multiplies each entry by a solid angle.
# i.e. (sum ∘ ⋆₀) computes a Riemann sum.
integral_of_curl(curl::VForm) = integral_of_curl(DualForm{2}(s0*curl.data))

# X is a primal 0, 
vort_ring(0.4, 6, PointVortexParams(3.0, 0.15), point_vortex)
u₀ = ComponentArray(dη = s0*X, β = 1e-9*s1*d0*X)

constants_and_parameters = (
  μ = 0.001,)
# TODO Units are probably heinous for u₀

@info "RMS of divergence of initial velocity: $(∘(RMS, div, curl_stream)(ψ))"
@info "Integral of initial curl: $(integral_of_curl(VForm(X)))"
end # ICs

###################

# η = dω::DualForm(1)
# DVF = map(sd[sd[:tri_center], :dual_point]) do (x,y,z); SVector(x, y^2, 0) end;

# # differential operator
# matt = ♭_mat(sd);

# # deRham map of our cohain
# primal1 = only.(matt*DVF);
# dual2 = dd1*(s1*primal1);

# β = map(sd[:point]) do (x,y,z); x^2 * 1e-3 end;

# u₀ = ComponentArray(dη = dual2, β = β, );

@info("Solving")
tₑ = 1.0; 

prob = ODEProblem(f, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob,
  Tsit5(),
  dtmax = 1e-3,
  dense=false,
  progress=true, progress_steps=1);

function save_gif(file_name, soln)
    time = Observable(0.0)
    fig = Figure()
    Label(fig[1, 1, Top()], @lift("...at $($time)"), padding = (0, 0, 5, 0))
    ax = CairoMakie.Axis(fig[1,1])
    msh = CairoMakie.mesh!(ax, s,
      color=@lift(s0inv*soln($time).dη),
      colormap=Reverse(:redsblues))
    Colorbar(fig[1,2], msh)
    record(fig, file_name, soln.t[1:10:end]; framerate = 10) do t
      time[] = t
    end
end
save_gif("vid3.mp4", soln)
 
