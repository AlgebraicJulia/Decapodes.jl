using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using CairoMakie
using ComponentArrays
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using StaticArrays
using Distributions

using LoggingExtras
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

lx = ly = 2π
dx = dy = 0.04 # 0.02
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
  map(x -> taylor_vortex(x, cntr, p), sd[:point])

function vort_ring(p::TaylorVortexParams, formula)
  sum(map(x -> formula(sd, x, p), [Point3d(lx / 2, ly / 2 - 0.4, 0.0), Point3d(lx / 2, ly / 2 + 0.4, 0.0)]))
end

pp♯ = ♯_mat(sd, AltPPSharp())
dd♯ = ♯_mat(sd, LLSDDSharp())

function plot_primal_vector_field(s, α)
  α♯ = pp♯ * α
  x_com = v -> dot(v, Point3d(1, 0, 0))
  y_com = v -> dot(v, Point3d(0, 1, 0))

  u = x_com.(α♯)
  v = y_com.(α♯)

  fig = Figure(size=(1000, 1000))
  ax_u = CairoMakie.Axis(fig[1, 1]; title="u", aspect=1)
  msh_u = CairoMakie.mesh!(ax_u, s, color=u, colormap=:jet)
  Colorbar(fig[1, 2], msh_u)

  ax_v = CairoMakie.Axis(fig[2, 1]; title="v", aspect=1)
  msh_v = CairoMakie.mesh!(ax_v, s, color=v, colormap=:jet)
  Colorbar(fig[2, 2], msh_v)

  colsize!(fig.layout, 1, Aspect(1, 1.0))
  resize_to_layout!(fig)

  return fig
end

function plot_dual_vector_field(s, α)
  α♯ = dd♯ * α
  x_com = v -> dot(v, Point3d(1, 0, 0))
  y_com = v -> dot(v, Point3d(0, 1, 0))

  u = interp * x_com.(α♯)
  v = interp * y_com.(α♯)

  fig = Figure(size=(1000, 1000))
  ax_u = CairoMakie.Axis(fig[1, 1]; title="u", aspect=1)
  msh_u = CairoMakie.mesh!(ax_u, s, color=u, colormap=:jet)
  Colorbar(fig[1, 2], msh_u)

  ax_v = CairoMakie.Axis(fig[2, 1]; title="v", aspect=1)
  msh_v = CairoMakie.mesh!(ax_v, s, color=v, colormap=:jet)
  Colorbar(fig[2, 2], msh_v)

  colsize!(fig.layout, 1, Aspect(1, 1.0))
  resize_to_layout!(fig)

  return fig
end

function plot_primal_vorticity(sd, α)
  α♯ = pp♯ * α

  filter_func = p -> (p[1] - pi)^2 + (p[2] - pi)^2 <= 1

  primal_p = sd[:point]
  primal_inner_points = findall(filter_func, primal_p)

  x = map(p -> p[1], primal_p[primal_inner_points])
  y = map(p -> p[2], primal_p[primal_inner_points])

  vel_u = map(p -> p[1], α♯[primal_inner_points])
  vel_v = map(p -> p[2], α♯[primal_inner_points])

  dual_p = sd[triangle_center(sd), :dual_point]
  dual_inner_points = findall(filter_func, dual_p)

  dual_x = map(p -> p[1], dual_p[dual_inner_points])
  dual_y = map(p -> p[2], dual_p[dual_inner_points])

  ω = abs.(hdg_2) * d1 * α

  fig = Figure()
  ax = CairoMakie.Axis(fig[1, 1]; title="Vorticity using the primal velocity")
  msh = CairoMakie.scatter!(ax, dual_x, dual_y, color=ω[dual_inner_points], colormap=:jet)
  msh = CairoMakie.arrows!(ax, x, y, vel_u, vel_v, lengthscale=0.1)
  display(fig)
end

function plot_dual_vorticity(sd, α)
  α♯ = dd♯ * α

  filter_func = p -> (p[1] - pi)^2 + (p[2] - pi)^2 <= 1

  primal_p = sd[:point]
  primal_inner_points = findall(filter_func, primal_p)

  x = map(p -> p[1], primal_p[primal_inner_points])
  y = map(p -> p[2], primal_p[primal_inner_points])

  dual_p = sd[triangle_center(sd), :dual_point]
  dual_inner_points = findall(filter_func, dual_p)

  dual_x = map(p -> p[1], dual_p[dual_inner_points])
  dual_y = map(p -> p[2], dual_p[dual_inner_points])

  dual_vel_u = map(p -> p[1], α♯[dual_inner_points])
  dual_vel_v = map(p -> p[2], α♯[dual_inner_points])

  ω = inv_hdg_0 * dual_d1 * α

  fig = Figure()
  ax = CairoMakie.Axis(fig[1, 1]; title="Vorticity using the dual velocity")
  msh = CairoMakie.scatter!(ax, x, y, color=ω[primal_inner_points], colormap=:jet)
  msh = CairoMakie.arrows!(ax, dual_x, dual_y, dual_vel_u, dual_vel_v, lengthscale=0.1)
  display(fig)
end

function plot_zeroform(s, f)
  isdual = length(f) .== ntriangles(s)

  fig = Figure()
  ax = CairoMakie.Axis(fig[1, 1])
  msh = CairoMakie.mesh!(ax, s, color=isdual ? interp * f : f, colormap=:jet)
  Colorbar(fig[1, 2], msh)
  display(fig)
end

# TODO: This interpolates the dual velocity field
function plot_velocity_diff(sd, dual_u, u)
  dual_u♯ = pp♯ * ♭♯ * dual_u
  primal_u♯ = pp♯ * u

  u♯ = primal_u♯ .- dual_u♯

  filter_func = p -> (p[1] - 3)^2 + (p[2] - 3)^2 <= 1

  primal_p = sd[:point]
  primal_inner_points = findall(filter_func, primal_p)

  x = map(p -> p[1], primal_p[primal_inner_points])
  y = map(p -> p[2], primal_p[primal_inner_points])

  vel_u = map(p -> p[1], u♯[primal_inner_points])
  vel_v = map(p -> p[2], u♯[primal_inner_points])

  fig = Figure()
  ax = CairoMakie.Axis(fig[1, 1]; title="Difference in velocity fields")
  CairoMakie.scatter!(ax, x, y; color=norm.(u♯))
  CairoMakie.arrows!(ax, x, y, vel_u, vel_v, lengthscale=1)
  display(fig)
end

ω = vort_ring(TaylorVortexParams(1.0, 0.3), taylor_vortex)
plot_zeroform(s, ω)

Δ0 = Δ(0, sd)
ψ = Δ0 \ ω
ψ .= ψ .- minimum(ψ)

plot_zeroform(s, ψ)

function p0_d0_interpolation(sd::HasDeltaSet2D)
  m = spzeros(ntriangles(sd), nv(sd))
  for tri in triangles(sd)
    tri_area = sd[tri, :area]
    for i in 0:5
      dual_tri = tri + i * ntriangles(sd)
      v, _ = dual_triangle_vertices(sd, dual_tri)
      dual_tri_area = sd[dual_tri, :dual_area]
      m[tri, v] += dual_tri_area / tri_area
    end
  end
  m
end

dual_interp = p0_d0_interpolation(sd)

dual_ψ = dual_interp * ψ
plot_zeroform(s, dual_ψ)

plot_zeroform(s, interp * dual_ψ .- ψ)

boundary_edges = findall(x -> x != 0, dual_d0 * ones(ntriangles(sd)))

dual_u = hdg_1 * d0 * ψ

# tmp = dual_d0 * dual_ψ
# tmp[boundary_edges] .= 0
# u = inv_hdg_1 * tmp

u = ♭♯ * dual_u

plot_primal_vector_field(s, u)
plot_dual_vector_field(s, dual_u)
plot_velocity_diff(sd, dual_u, u)

# TODO: Primal is incorrect while the dual is right
plot_primal_vorticity(sd, u)
plot_dual_vorticity(sd, dual_u)

plot_zeroform(s, abs.(hdg_2) * d1 * u)

κ = 280 * 300 # R (dry gas constant) * T = 300K, P=ρRVT
ρ₀ = ones(nv(sd))
U₀ = wdg_10(u, ρ₀)

#                  d0                   d1
#       P0         ->     P1           ->     P2
#  hdg0 | inv_hdg_0  hdg1 | inv_hdg_1   hdg2 | inv_hdg_2
#       D2         <-     D1       <-         D0
#                dual_d1        dual_d0
#
# (U,u)::Form1
# (ρ, P)::Form0

form_one_interp = ♭♯ * hdg_1
form_zero_interp = interp * abs.(hdg_2) * d1
form_two_interp = ♭♯ * dual_d0 * abs.(hdg_2)

lap_first_term = inv_hdg_1 * dual_d0 * abs.(hdg_2) * d1
lap_second_term = d0 * inv_hdg_0 * dual_d1 * hdg_1

function momentum_continuity(U, ρ)
  u = wdg_10(U, 1 ./ ρ)

  return -wdg_10(U, codif_1 * u) - # U ∧ δu
         form_two_interp * wdg_11(u, form_one_interp * U) - # L(u, U)
         form_one_interp * wdg_10(u, form_zero_interp * U) +
         0.5 * wdg_10(form_two_interp * wdg_11(u, form_one_interp * u), ρ) + # 1/2 * ρ * d||u||^2
         d0 * (κ * ρ) + # dP, P = κρ
         μ * (lap_first_term * u + lap_second_term * u)# μΔu

end

function run_compressible_ns(U₀, ρ₀, tₑ, Δt)

  U = deepcopy(U₀)
  ρ = deepcopy(ρ₀)

  U_half = zeros(ne(sd))

  ρ_half = zeros(nv(sd))
  ρ_full = zeros(nv(sd))


  steps = ceil(Int64, tₑ / Δt)

  Us = [deepcopy(U₀)]
  rhos = [deepcopy(ρ₀)]

  for step in 1:steps
    U_half .= U .+ Δt / 2 * momentum_continuity(U, ρ)
    ρ_full .= ρ + Δt * codif_1 * U_half
    ρ_half .= 0.5 * (ρ + ρ_full)

    U .= U .+ Δt * momentum_continuity(U_half, ρ_half)
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

    push!(Us, deepcopy(U))
    push!(rhos, deepcopy(ρ))

    if step % 1000 == 0
      println("Loading simulation results: $(step / steps * 100)%")
    end
  end

  return Us, rhos
end

μ = 0 # \Re

tₑ = 0.5
Δt = dx / 360
Us, rhos = run_compressible_ns(U₀, ρ₀, tₑ, Δt);

ω_end = abs.(hdg_2) * d1 * wdg_10(Us[end], 1 ./ rhos[end])

dp = sd[triangle_center(sd), :dual_point]
x = map(p -> p[1], dp)
y = map(p -> p[2], dp)

scatter(x, y; color=ω_end)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ω_end, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

u₀ = wdg_10(Us[end], 1 ./ rhos[end])

plot_primal_vorticity(sd, u₀)

ω_init = hdg_2 * d1 * wdg_10(U₀, 1 ./ ρ₀)

dp = sd[triangle_center(sd), :dual_point]
x = map(p -> p[1], dp)
y = map(p -> p[2], dp)

scatter(x, y; color=ω_init, colormap=:jet)

ω_init = interp * hdg_2 * d1 * wdg_10(U₀, 1 ./ ρ₀)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ω_init, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=rhos[end], colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

