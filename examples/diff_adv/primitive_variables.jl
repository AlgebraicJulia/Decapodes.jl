using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using CairoMakie
using ComponentArrays
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays

using LoggingExtras
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

lx = ly = 2π
dx = dy = 0.08
s = triangulated_grid(lx, ly, dx, dy, Point3d, false)
nx = length(0:dx:lx)
ny = length(0:dy:ly)

sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

d0 = dec_differential(0, sd)
d1 = dec_differential(1, sd)

dual_d0 = dec_dual_derivative(0, sd);
dual_d1 = dec_dual_derivative(1, sd);

dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));

hdg_1 = dec_hodge_star(1, sd, DiagonalHodge())
hdg_2 = dec_hodge_star(2, sd, DiagonalHodge())

inv_hdg_0 = dec_inv_hodge_star(0, sd, DiagonalHodge())
inv_hdg_1 = dec_inv_hodge_star(1, sd, DiagonalHodge())

wdg_01 = dec_wedge_product(Tuple{0,1}, sd)

♭♯ = ♭♯_mat(sd)

div = hdg_2 * d1 * inv_hdg_1
dΔ0 = div * dual_d0
fdΔ0 = factorize(dΔ0);

Wv(v) = 0.5 * spdiagm(v) * abs.(d0)

struct TaylorVortexParams
  G::Real
  a::Real
end

function taylor_vortex(pnt::Point3d, cntr::Point3d, p::TaylorVortexParams)
  r = norm(pnt .- cntr)
  (p.G / p.a) * (2 - (r / p.a)^2) * exp(0.5 * (1 - (r / p.a)^2))
end

taylor_vortex(sd::HasDeltaSet, cntr::Point3d, p::TaylorVortexParams) =
  map(x -> taylor_vortex(x, cntr, p), point(sd))

function vort_ring(p::TaylorVortexParams, formula)
  sum(map(x -> formula(sd, x, p), [Point3d(lx / 2 - 0.4, ly / 2, 0.0), Point3d(lx / 2 + 0.4, ly / 2, 0.0)]))
end

ω = vort_ring(TaylorVortexParams(1.0, 0.3), taylor_vortex)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ω, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

Δ0 = Δ(0, sd)

ψ = Δ0 \ ω
ψ .= ψ .- minimum(ψ)
# u = hdg_1 * d0 * ψ
u_star = d0 * ψ

Δt = 1e-3
tₑ = 10
μ = 0

nt = ntriangles(sd)

rhs_U_mat = (-1 / Δt) * I + μ * d0 * inv_hdg_0 * dual_d1 * hdg_1
rhs_Pd_mat = inv_hdg_1 * dual_d0

function generate_F(u_star)
  v = ♭♯ * hdg_1 * u_star
  return (-1 / Δt) * u_star - μ * d0 * inv_hdg_0 * dᵦ * v + # Diffusion of v
         Wv(v) * inv_hdg_0 * (dual_d1 * hdg_1 * u_star + dᵦ * v) # Advection 
end

F₁ = generate_F(u_star)
F₂ = zeros(nt)

F = vcat(F₁, F₂)

rhs_top = hcat(rhs_U_mat, rhs_Pd_mat)
rhs_bottom = hcat(d1, zeros(nt, nt))

rhs = vcat(rhs_top, rhs_bottom)

f_rhs = factorize(rhs)

U = zeros(ne(sd) + nt)
steps = ceil(Int64, tₑ / Δt)
for step in 1:steps
  F = vcat(generate_F(u_star), F₂)
  U .= f_rhs \ F
  u_star .= U[1:ne(sd)]

  if step % 1000 == 0
    println("Loading simulation results: $(step / steps * 100)%")
  end
end

ω_end = inv_hdg_0 * dual_d1 * hdg_1 * u_star

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ω_end, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

# u₀ = ComponentArray(u=u)

# mutable struct Solution
#   t::Vector{Float64}
#   u::Vector{ComponentVector{Float64}}

#   function Solution(t₀::Float64, u₀::ComponentVector{Float64})
#     return new(Float64[t₀], ComponentVector{Float64}[u₀])
#   end
# end

# function euler_equation_with_projection(u₀, tₑ, Δt)

#   soln = Solution(0.0, u₀)
#   u = deepcopy(u₀)

#   u_int = zeros(ne(sd))
#   p_next = zeros(ntriangles(sd))

#   steps = ceil(Int64, tₑ / Δt)
#   for step in 1:steps
#     v = ♭♯ * u
#     u_int .= u .- hdg_1 * wdg_01(inv_hdg_0 * dual_d1 * u, v) * Δt
#     p_next .= fdΔ0 \ (div * u_int) / Δt
#     u .= u_int .- Δt * (dual_d0 * p_next)

#     push!(soln.t, 0 + step * Δt)
#     push!(soln.u, ComponentVector(u=u))
#   end
#   return soln
# end

# soln = euler_equation_with_projection(u₀, tₑ, Δt)

# ω_end = inv_hdg_0 * dual_d1 * soln.u[end]

# fig = Figure();
# ax = CairoMakie.Axis(fig[1, 1])
# msh = CairoMakie.mesh!(ax, s, color=ω_end, colormap=:jet)
# Colorbar(fig[1, 2], msh)
# display(fig)

