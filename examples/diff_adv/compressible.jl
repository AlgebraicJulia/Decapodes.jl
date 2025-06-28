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
dx = dy = 0.08 # 0.02
s = triangulated_grid(lx, ly, dx, dy, Point3d, false)
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
inv_hdg_2 = dec_inv_hodge_star(2, sd, DiagonalHodge())

wdg_10 = dec_wedge_product(Tuple{1,0}, sd)
wdg_11 = dec_wedge_product(Tuple{1,1}, sd)

interp = d0_p0_interpolation(sd; hodge=DiagonalHodge())
♭♯ = ♭♯_mat(sd)

codif_1 = inv_hdg_0 * dual_d1 * hdg_1

struct TaylorVortexParams
  G::Real
  a::Real
end

function taylor_vortex(pnt::Point3d, cntr::Point3d, p::TaylorVortexParams)
  r = norm(pnt .- cntr)
  (p.G / p.a) * (2 - (r / p.a)^2) * exp(0.5 * (1 - (r / p.a)^2))
end

taylor_vortex(sd::HasDeltaSet, cntr::Point3d, p::TaylorVortexParams) =
  map(x -> taylor_vortex(x, cntr, p), sd[triangle_center(sd), :dual_point])

function vort_ring(p::TaylorVortexParams, formula)
  sum(map(x -> formula(sd, x, p), [Point3d(lx / 2 - 0.4, ly / 2, 0.0), Point3d(lx / 2 + 0.4, ly / 2, 0.0)]))
end

ω = vort_ring(TaylorVortexParams(1.0, 0.3), taylor_vortex)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=interp * ω, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

# Δ0 = Δ(0, sd)
# ψ = Δ0 \ ω
# ψ .= ψ .- minimum(ψ)
# u = ♭♯ * hdg_1 * d0 * ψ
# u = hdg_1 * d0 * ψ

dΔ0 = hdg_2 * d1 * inv_hdg_1 * dual_d0
fdΔ0 = factorize(dΔ0);
ψ = dΔ0 \ ω
ψ .= ψ .- minimum(ψ)
u = inv_hdg_1 * dual_d0 * ψ

κ = 280 * 300 # R (dry gas constant) * T = 300K, P=ρRVT
ρ₀ = ones(nv(sd))
U₀ = wdg_10(u, ρ₀)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ρ₀, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

#                  d0                   d1
#       P0         ->     P1           ->     P2
#  hdg0 | inv_hdg_0  hdg1 | inv_hdg_1   hdg2 | inv_hdg_2
#       D2         <-     D1       <-         D0
#                dual_d1        dual_d0
#
# (U,u)::Form1
# (ρ, P)::Form0

form_one_interp = ♭♯ * hdg_1
form_zero_interp = interp * hdg_2 * d1
form_two_interp = ♭♯ * dual_d0 * hdg_2

lap_first_term = inv_hdg_1 * dual_d0 * hdg_2 * d1
lap_second_term = d0 * inv_hdg_0 * dual_d1 * hdg_1

function momentum_continuity(U, ρ, Δt)
  u = wdg_10(U, 1 ./ ρ)

  return -wdg_10(U, codif_1 * u) - # U ∧ δu
         form_two_interp * wdg_11(u, form_one_interp * U) - # L(u, U)
         form_one_interp * wdg_10(u, form_zero_interp * U) +
         0.5 * wdg_10(form_two_interp * wdg_11(u, form_one_interp * u), ρ) + # 1/2 * ρ * d||u||^2
         d0 * (κ * ρ) + # dP, P = κρ
         μ * (lap_first_term * u + lap_second_term * u) # μΔu

end

function run_compressible_ns(U₀, ρ₀, tₑ, Δt)

  U = deepcopy(U₀)
  ρ = deepcopy(ρ₀)

  U_half = zeros(ne(sd))
  U_full = zeros(ne(sd))

  ρ_half = zeros(nv(sd))
  ρ_full = zeros(nv(sd))


  steps = ceil(Int64, tₑ / Δt)

  Us = [U₀]
  rhos = [ρ₀]

  for step in 1:steps
    U_half .= U .+ Δt / 2 * momentum_continuity(U, ρ, Δt)
    ρ_full .= ρ + Δt * codif_1 * U_half
    ρ_half .= 0.5 * (ρ + ρ_full)

    U .= U .+ Δt * momentum_continuity(U_half, ρ_half, Δt)
    ρ .= ρ_full

    if any(isnan.(U))
      println("Warning, NAN result in U at step: $(step)")
      break
    elseif any(isinf.(U))
      println("Warning, INF result in U at step: $(step)")
      break
    elseif any(isnan.(ρ))
      println("Warning, NAN result in ρ at step: $(step)")
      break
    elseif any(isinf.(ρ))
      println("Warning, INF result in ρ at step: $(step)")
      break
    end

    push!(Us, U)
    push!(rhos, ρ)

    if step % 1000 == 0
      println("Loading simulation results: $(step / steps * 100)%")
    end
  end

  return Us, rhos
end

μ = 0 # \Re

tₑ = 1
Δt = 1e-4
Us, rhos = run_compressible_ns(U₀, ρ₀, tₑ, Δt)

timestep = 8
ω_end = interp * hdg_2 * d1 * wdg_10(Us[end], 1 ./ rhos[end])

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ω_end, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=rhos[end], colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)
