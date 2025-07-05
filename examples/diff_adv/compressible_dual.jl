using CairoMakie
using CombinatorialSpaces
using ComponentArrays
using LinearAlgebra
using MLStyle
using SparseArrays
using StaticArrays

import CombinatorialSpaces.FastDEC: wedge_dd_01_mat

using CUDA
using CUDA.CUSPARSE

lx = ly = 2π
dx = dy = 0.08
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

wdg_10 = dec_wedge_product(Tuple{1,0}, sd, Val{:CUDA}, CuArray, Float64)
wdg_11 = dec_wedge_product(Tuple{1,1}, sd, Val{:CUDA}, CuArray, Float64)

# TODO: Running d0_p0 directly is causing OutOfMemory errors
# interp = d0_p0_interpolation(sd; hodge=DiagonalHodge())
interp = SparseMatrixCSC{Float64}(inv_hdg_0) * p2_d2_interpolation(sd) * SparseMatrixCSC{Float64}(inv_hdg_2)
♭♯ = ♭♯_mat(sd)

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
  sum(map(x -> formula(sd, x, p), [Point3d(lx / 2 - 0.4, ly / 2, 0.0), Point3d(lx / 2 + 0.4, ly / 2, 0.0)]))
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

boundary_edges = findall(x -> x != 0, dual_d0 * ones(ntriangles(sd)))

u = hdg_1 * d0 * ψ

plot_dual_vector_field(s, u)
plot_dual_vorticity(sd, u)

plot_zeroform(s, inv_hdg_0 * dual_d1 * u)

wdg_10_dd_mat = CuSparseMatrixCSC{Float64}(wedge_dd_01_mat(sd))

wdg_10dd(α, f) = α .* (wdg_10_dd_mat * f) 

κ = 280 * 300 # R (dry gas constant) * T = 300K, P=ρRVT
ρ₀ = CUDA.ones(Float64, ntriangles(sd))
U₀ = wdg_10dd(CuArray{Float64}(u), ρ₀)

#                  d0                   d1
#       P0         ->     P1           ->     P2
#  hdg0 | inv_hdg_0  hdg1 | inv_hdg_1   hdg2 | inv_hdg_2
#       D2         <-     D1       <-         D0
#                dual_d1        dual_d0
#
# (U,u)::DualForm1
# (ρ, P)::DualForm0

m₀ = CUDA.sum(ρ₀)

form_one_interp_1 = CuSparseMatrixCSC{Float64}(inv_hdg_1)
form_one_interp_2 = CuSparseMatrixCSC{Float64}(hdg_1)
form_zero_interp = CuSparseMatrixCSC{Float64}(inv_hdg_0 * dual_d1)
form_two_interp = CuSparseMatrixCSC{Float64}(dual_d0 * hdg_2)

special_dual_d0 = deepcopy(dual_d0)
special_dual_d0[boundary_edges, :] .= 0
cu_dual_d0 = CuSparseMatrixCSC{Float64}(special_dual_d0)

codif_1 = CuSparseMatrixCSC{Float64}(abs.(hdg_2)) * CuSparseMatrixCSC{Float64}(d1) * CuSparseMatrixCSC{Float64}(inv_hdg_1)

lap_first_term = CuSparseMatrixCSC{Float64}(hdg_1 * d0 * inv_hdg_0 * dual_d1)
lap_second_term = CuSparseMatrixCSC{Float64}(special_dual_d0 * hdg_2 * d1 * inv_hdg_1)

lap_term = lap_first_term + lap_second_term

function momentum_continuity(U, ρ)
  u = wdg_10dd(U, 1 ./ ρ)

  return -wdg_10dd(U, codif_1 * u) - # U ∧ δu
         form_two_interp * wdg_11(u, form_one_interp_1 * U) - # L(u, U)
         form_one_interp_2 * wdg_10(u, form_zero_interp * U) +
         0.5 * wdg_10dd(form_two_interp * wdg_11(u, form_one_interp_1 * u), ρ) + # 1/2 * ρ * d||u||^2
         cu_dual_d0 * (κ * ρ) + # dP, P = κρ
         μ * lap_term * u # μΔu

end

function run_compressible_ns(U₀, ρ₀, tₑ, Δt; saveat=500)

  U = deepcopy(U₀)
  ρ = deepcopy(ρ₀)

  U_half = CUDA.zeros(Float64, ne(sd))

  ρ_half = CUDA.zeros(Float64, ntriangles(sd))
  ρ_full = CUDA.zeros(Float64, ntriangles(sd))


  steps = ceil(Int64, tₑ / Δt)

  Us = [deepcopy(U₀)]
  rhos = [deepcopy(ρ₀)]

  for step in 1:steps
    U_half .= U .+ 0.5 * Δt * momentum_continuity(U, ρ)
    ρ_full .= ρ + Δt * codif_1 * U_half
    ρ_half .= 0.5 .* (ρ .+ ρ_full)

    U .= U .+ Δt * momentum_continuity(U_half, ρ_half)
    ρ .= ρ_full

    if any(CUDA.isnan.(U))
      println("Warning, NAN result in U at step: $(step)")
      break
    elseif any(CUDA.isinf.(U))
      println("Warning, INF result in U at step: $(step)")
      break
    elseif any(CUDA.isnan.(ρ))
      println("Warning, NAN result in ρ at step: $(step)")
      break
    elseif any(CUDA.isinf.(ρ))
      println("Warning, INF result in ρ at step: $(step)")
      break
    end

    if step % saveat == 0
      push!(Us, deepcopy(U))
      push!(rhos, deepcopy(ρ))
      println("Loading simulation results: $((step / steps) * 100)%")
      println("Relative mass is : $((CUDA.sum(ρ) / m₀) * 100)%")
      println("-----")
    end
  end

  return Us, rhos
end

μ = 0 # \Re

tₑ = 1
Δt = 1e-4 # dx / 360
Us, rhos = run_compressible_ns(U₀, ρ₀, tₑ, Δt; saveat=500);

timestep = length(Us)
u_end = Array(wdg_10dd(Us[timestep], 1 ./ rhos[timestep]))
ω_end = inv_hdg_0 * dual_d1 * u_end

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=ω_end, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

# plot_primal_vorticity(sd, u_end)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s, color=interp*Array(rhos[timestep]), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

function save_dynamics(save_file_name)
  time = Observable(1)

  w = @lift(inv_hdg_0 * dual_d1 * Array(wdg_10dd(Us[$time], 1 ./ rhos[$time])))
  ρ = @lift(interp*Array(rhos[$time]))

  f = Figure(size=(1000, 1000))

  ax_w = CairoMakie.Axis(f[1, 1], title="Vorticity with μ=$(μ)")
  msh_w = mesh!(ax_w, s; color=w, colormap=:jet, colorrange=extrema(ω))
  Colorbar(f[1, 2], msh_w)

  ax_ρ = CairoMakie.Axis(f[2, 1], title="Density with μ=$(μ)")
  msh_ρ = mesh!(ax_ρ, s; color=ρ, colormap=:jet, colorrange=(-1e-6,1e-6).+1)
  Colorbar(f[2, 2], msh_ρ)

  colsize!(fig.layout, 1, Aspect(1, 1.0))
  resize_to_layout!(fig)

  timestamps = 1:length(Us)
  CairoMakie.record(f, save_file_name, timestamps; framerate=15) do t
    time[] = t
  end
end

save_dynamics("Compressible_Taylor_Vortices.mp4")
