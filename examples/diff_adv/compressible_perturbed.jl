using CairoMakie
using CombinatorialSpaces
using ComponentArrays
using Distributions
using LinearAlgebra
using MLStyle
using SparseArrays
using StaticArrays
using Printf
using JLD2

using CUDA
using CUDA.CUSPARSE

# lx = ly = 8
lx = 4; ly = 8
dx = dy = 0.02
s = triangulated_grid(lx, ly, dx, dy, Point3d, false)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

d0 = dec_differential(0, sd);
d1 = dec_differential(1, sd);

dual_d0 = dec_dual_derivative(0, sd);
dual_d1 = dec_dual_derivative(1, sd);

hdg_1 = dec_hodge_star(1, sd, DiagonalHodge());
hdg_2 = dec_hodge_star(2, sd, DiagonalHodge());

inv_hdg_0 = dec_inv_hodge_star(0, sd, DiagonalHodge());
inv_hdg_1 = dec_inv_hodge_star(1, sd, DiagonalHodge());
inv_hdg_2 = dec_inv_hodge_star(2, sd, DiagonalHodge());

codif_1 = SparseMatrixCSC(inv_hdg_0) * dual_d1 * SparseMatrixCSC(hdg_1);

wdg_10 = dec_wedge_product(Tuple{1,0}, sd, Val{:CUDA}, CuArray, Float64);
wdg_11 = dec_wedge_product(Tuple{1,1}, sd, Val{:CUDA}, CuArray, Float64);

# TODO: Running d0_p0 directly is causing OutOfMemory errors
interp = d0_p0_interpolation(sd; hodge=DiagonalHodge())
flatsharp = ♭♯_mat(sd);

pp♭ = ♭_mat(sd, PPFlat());

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

# dual_interp = p0_d0_interpolation(sd)

pp♯ = ♯_mat(sd, AltPPSharp());
dd♯ = ♯_mat(sd, LLSDDSharp());

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

function plot_zeroform(s, f; title=nothing, colormap=:jet)
  isdual = length(f) .== ntriangles(s)

  fig = Figure()
  ax = isnothing(title) ? CairoMakie.Axis(fig[1, 1]) : CairoMakie.Axis(fig[1, 1]; title=title)
  msh = CairoMakie.mesh!(ax, s, color=isdual ? interp * f : f, colormap=colormap)
  Colorbar(fig[1, 2], msh)
  display(fig)
  return fig
end

# TODO: This interpolates the dual velocity field
function plot_velocity_diff(sd, dual_u, u)
  dual_u♯ = pp♯ * flatsharp * dual_u
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

abstract type AbstractBoundaryMapping end

struct TriangleMapping <: AbstractBoundaryMapping
  boundary_zone::AbstractVector{Int64}
  copy_to_boundary_zone::AbstractVector{Int64}
end

struct VertexMapping <: AbstractBoundaryMapping
  boundary_zone::AbstractVector{Int64}
  copy_to_boundary_zone::AbstractVector{Int64}

  function VertexMapping(bz::AbstractVector{Int64}, cbz::AbstractVector{Int64})
    new(bz, cbz)
  end
end

struct EdgeMapping <: AbstractBoundaryMapping
  boundary_zone::AbstractVector{Int64}
  copy_to_boundary_zone::AbstractVector{Int64}
end

function VertexMapping(sd, t::TriangleMapping)
  vertex_boundary = vcat(triangle_vertices(sd, t.boundary_zone)...)
  vertex_copy_to = vcat(triangle_vertices(sd, t.copy_to_boundary_zone)...)
  return VertexMapping(vertex_boundary, vertex_copy_to)
end

function EdgeMapping(sd, t::TriangleMapping)
  edge_boundary = vcat(triangle_edges(sd, t.boundary_zone)...)
  edge_copy_to = vcat(triangle_edges(sd, t.copy_to_boundary_zone)...)
  return EdgeMapping(edge_boundary, edge_copy_to)
end

function collect_bottom_boundary_triangles(dx, lx, dy, ly, depth)
  nx = length(0:dx:lx)
  ny = length(0:dy:ly)

  boundary_triangles = Int64[]
  row = 2(nx - 1)
  for y in 1:2*depth
    append!(boundary_triangles, collect(1:2(nx-1)) .+ (y - 1) * row)
  end
  return TriangleMapping(boundary_triangles, boundary_triangles .+ (ny - 4 * depth - 1) * row)
end

function collect_top_boundary_triangles(dx, lx, dy, ly, depth)
  nx = length(0:dx:lx)
  ny = length(0:dy:ly)

  boundary_triangles = Int64[]
  row = 2(nx - 1)
  for y in 1:2*depth
    append!(boundary_triangles, collect(1:2(nx-1)) .+ (ny - y - 1) * row)
  end
  return TriangleMapping(boundary_triangles, boundary_triangles .- (ny - 4 * depth - 1) * row)
end

function collect_left_boundary_triangles(dx, lx, dy, ly, depth)
  nx = length(0:dx:lx)
  ny = length(0:dy:ly)

  boundary_triangles = Int64[]
  row = 2(nx - 1)
  xs = collect(1:2depth)
  for y in 1:ny-1
    append!(boundary_triangles, xs .+ (y - 1) * row)
  end
  return TriangleMapping(boundary_triangles, boundary_triangles .+ (row .- 4 .* depth))
end

function collect_right_boundary_triangles(dx, lx, dy, ly, depth)
  nx = length(0:dx:lx)
  ny = length(0:dy:ly)

  boundary_triangles = Int64[]
  row = 2(nx - 1)
  xs = collect(1:2depth)
  for y in 1:ny-1
    append!(boundary_triangles, y * row .- xs .+ 1)
  end
  return TriangleMapping(boundary_triangles, boundary_triangles .- (row .- 4 .* depth))
end

depth = 5
topb = collect_top_boundary_triangles(dx, lx, dy, ly, depth);
botb = collect_bottom_boundary_triangles(dx, lx, dy, ly, depth);

leftb = collect_left_boundary_triangles(dx, lx, dy, ly, depth);
rightb = collect_right_boundary_triangles(dx, lx, dy, ly, depth);

import Base.vcat
vcat(b1::AbstractBoundaryMapping, b2::AbstractBoundaryMapping) = TriangleMapping(vcat(b1.boundary_zone, b2.boundary_zone),
  vcat(b1.copy_to_boundary_zone, b2.copy_to_boundary_zone))

function apply_periodic!(x, m::AbstractBoundaryMapping)
  x[m.boundary_zone] .= x[m.copy_to_boundary_zone]
  return x
end

tb_bounds = vcat(topb, botb);
lr_bounds = vcat(leftb, rightb);
vtb_bounds = VertexMapping(sd, tb_bounds);
vlr_bounds = VertexMapping(sd, lr_bounds);
etb_bounds = EdgeMapping(sd, tb_bounds);
elr_bounds = EdgeMapping(sd, lr_bounds);

##############################
##### INITIAL CONDITIONS #####
##############################

const P₀ = 1e5 # Reference pressure for potential temperature
const Cₚ = 1006 # Specific heat
const R = 287 # Specific gas constant
const ρᵣ = 1 # Reference density
const θᵣ = 300 # Reference potential temperature, temperature at P₀
const Pᵣ = ρᵣ * R * θᵣ
const gₐ = -9.81 # Acceleration due to gravity
g = CuVector{Float64}(eval_constant_primal_form(sd, Point3d(0, gₐ, 0)))

# TODO: This assumes theta is constant in z
function hydrostatic_pressure(theta, h)
  return (1 ./ (Cₚ .* theta) .* (P₀)^(R / Cₚ) .* gₐ .* h .+ Pᵣ^(R / Cₚ)) .^ (Cₚ / R)
end

function hydrostatic_density(theta, h)
  Pₕ = hydrostatic_pressure(theta, h)

  return (Pₕ / (R * theta)) * (P₀ / Pₕ)^(R / Cₚ)
end

hydrostatic_pressure(300, 0) ≈ Pᵣ
hydrostatic_density(300, 0)

function generate_hydrostatic_vars(sd, theta=300)
  thetaₕ = theta * CUDA.ones(Float64, nv(sd))
  ρₕ = CuVector{Float64}(map(p -> hydrostatic_density(theta, p[2]), sd[:point]))
  return ρₕ .* thetaₕ, ρₕ
end

# From the ideal gas law
function exact_pressure(Theta)
  R_Cₚ = R / Cₚ
  return (Theta .* R .* (P₀ .^ -R_Cₚ)) .^ (1 / (1 - R_Cₚ))
end

# XXX: Expect reasonable results with +-30 of hydrostatic
function perturbed_pressure(Thetaₕ, Theta′)
  R_Cₚ = R / Cₚ
  r_R_Cₚ = 1 / (1 - R_Cₚ)
  return r_R_Cₚ .* (Thetaₕ) .^ (r_R_Cₚ - 1) .* Theta′ .* (R .* (P₀ .^ -R_Cₚ)) .^ (r_R_Cₚ)
end

# # HYDROSTATIC
# simname = "Perturbed_Hydrostatic"
# Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)
# U₀ = CUDA.zeros(Float64, ne(sd))
# ρ′ = CUDA.zeros(Float64, nv(sd))
# Theta′ = CUDA.zeros(Float64, nv(sd))

# # WARM BUBBLE
# simname = "Perturbed_Warm_Bubble"
# Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)
# U₀ = CUDA.zeros(Float64, ne(sd))

# ρ′ = CUDA.zeros(Float64, nv(sd))
# theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
# Theta′ = ρₕ .* CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])

# # WARM BUBBLE
# simname = "Perturbed_Warm_Bubble_Gradient"
# thetaₕ = 300 * CUDA.ones(Float64, nv(sd))
# Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)
# Pₕ = exact_pressure(Thetaₕ)

# perturb = false
# peturbation_size = 1e-4

# U₀ = CUDA.zeros(Float64, ne(sd))
# theta_dist = MvNormal([lx / 2, ly / 4], 0.25)
# theta′ = 10 * CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]]) .+ perturb .* peturbation_size .* (2 * CUDA.rand(Float64, nv(sd)) .- 1)
# plot_zeroform(s, Array(theta′))

# ρ′ = (Pₕ .^ (1 - R / Cₚ) * P₀^(R / Cₚ)) ./ (R .* (thetaₕ .+ theta′)) .- ρₕ .+ perturb .* peturbation_size .* (2 * CUDA.rand(Float64, nv(sd)) .- 1)
# plot_zeroform(s, Array(ρ′))
# Theta′ = ρ′ .* thetaₕ .+ ρₕ .* theta′ .+ ρ′ .* theta′
# plot_zeroform(s, Array(Theta′))

# # WARM BUBBLE WITH WIND
# simname = "Perturbed_Warm_Bubble_Wind"
# Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)

# ρ′ = CUDA.zeros(Float64, nv(sd))
# theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
# Theta′ = ρₕ .* CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])

# U₀ = (ρₕ .+ ρ′) .* CuVector{Float64}(eval_constant_primal_form(sd, Point3d(1, 0, 0)))

# # COLD BUBBLE
# simname = "Perturbed_Cold_Bubble"
# Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)
# U₀ = CUDA.zeros(Float64, ne(sd))

# ρ′ = CUDA.zeros(Float64, nv(sd))
# theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
# Theta′ = ρₕ .* -CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])

# # LIGHT BUBBLE
# simname = "Perturbed_Light_Bubble_Gradient"
# thetaₕ = 300 * CUDA.ones(Float64, nv(sd))
# Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)
# Pₕ = exact_pressure(Thetaₕ)

# U₀ = CUDA.zeros(Float64, ne(sd))
# rho_dist = MvNormal([lx / 2, ly / 4], 0.5)
# ρ′ = -0.01 * CuArray{Float64}([pdf(rho_dist, [p[1], p[2]]) for p in sd[:point]])

# theta′ = (Pₕ .^ (1 - R / Cₚ) * P₀^(R / Cₚ)) ./ (R .* (ρₕ .+ ρ′)) .- 300
# plot_zeroform(s, Array(theta′))
# Theta′ = ρ′ .* thetaₕ .+ ρₕ .* theta′ .+ ρ′ .* theta′
# plot_zeroform(s, Array(Theta′))

# RAYLEIGH-TAYLOR INSTABILITY
simname = "Perturbed_RTI"
thetaₕ = 300 * CUDA.ones(Float64, nv(sd))
Thetaₕ, ρₕ = generate_hydrostatic_vars(sd)
Pₕ = exact_pressure(Thetaₕ)

perturb = true
peturbation_size = 1e-3

k = 4 / lx
U_dist = MvNormal([lx / 2, ly / 2], [0.5, 0.5])
U₀ = CuArray{Float64}(only.(pp♭ * map(p -> SVector{3, Float64}(0, 2 * pdf(U_dist, [p[1], p[2]]) * sin(2pi*k*p[1]), 0), sd[:point]))) # CUDA.zeros(Float64, ne(sd))
C = 10; D = -20
theta′ = CuArray{Float64}(map(p -> D*tanh(C*(p[2] - ly/2)), sd[:point])) .+ perturb .* peturbation_size .* (2 * CUDA.rand(Float64, nv(sd)) .- 1)
plot_zeroform(s, Array(theta′))

ρ′ = (Pₕ .^ (1 - R / Cₚ) * P₀^(R / Cₚ)) ./ (R .* (thetaₕ .+ theta′)) .- ρₕ .+ perturb .* peturbation_size .* (2 * CUDA.rand(Float64, nv(sd)) .- 1)
plot_zeroform(s, Array(ρ′))
Theta′ = ρ′ .* thetaₕ .+ ρₕ .* theta′ .+ ρ′ .* theta′
plot_zeroform(s, Array(Theta′))

plot_primal_vector_field(s, Array(U₀))

##################################
##### END INITIAL CONDITIONS #####
##################################

plot_zeroform(s, Array(ρₕ))
plot_zeroform(s, Array(Thetaₕ))
plot_zeroform(s, Array(exact_pressure(Thetaₕ)))

plot_zeroform(s, Array(ρ′))
plot_zeroform(s, Array(Theta′))
plot_zeroform(s, Array(perturbed_pressure(Thetaₕ, Theta′)))

# function smoothing(sd, c_smooth)
#   mat = spzeros(nv(sd), nv(sd))
#   for e in edges(sd)
#     v1 = sd[e, :∂v0]
#     v2 = sd[e, :∂v1]
#     w = 1 / sd[e, :length]
#     mat[v1, v2] = w
#     mat[v2, v1] = w
#   end

#   c = c_smooth ./ 2

#   for v in vertices(sd)
#     row = mat[v, :]
#     tot_w = sum(row)
#     for i in row.nzind
#       mat[v, i] = c .* row[i] ./ tot_w
#     end

#     mat[v, v] = (1 - c)
#   end
#   return mat
# end

function smoothing(sd, c_smooth)
  n = nv(sd) + 2*ne(sd)
  I = zeros(Int32, n)
  J = zeros(Int32, n)
  V = zeros(Float64, n) # Positive smoothing
  V2 = zeros(Float64, n) # Negative smoothing

  c = c_smooth / 2

  aggr = zeros(nv(sd))
  for e in edges(sd)
    i = 2*e - 1

    v1 = sd[e, :∂v0]
    v2 = sd[e, :∂v1]
    w = 1 / sd[e, :length]

    I[i] = v1; J[i] = v2; V[i] = w
    I[i+1] = v2; J[i+1] = v1; V[i+1] = w

    aggr[v1] += w
    aggr[v2] += w
  end

  for i in 1:2ne(sd)
    V[i] = c * V[i] / aggr[I[i]]
    V2[i] = -V[i]
  end

  for v in vertices(sd)
    i = 2*ne(sd) + v
    I[i] = v; J[i] = v; V[i] = 1 - c; V2[i] = 1 + c
  end

  return sparse(I, J, V2) * sparse(I, J, V)
end


Theta_smooth = 0.2
Theta_smoothing_mat = CuSparseMatrixCSC{Float64}(smoothing(sd, Theta_smooth));

rho_smooth = 0.05
rho_smoothing_mat = CuSparseMatrixCSC{Float64}(smoothing(sd, rho_smooth));

momentum_flatsharp_hdg_1 = CuSparseMatrixCSC{Float64}(flatsharp * hdg_1);
momentum_interp_hdg2_d1 = CuSparseMatrixCSC{Float64}(interp * hdg_2 * d1);
momentum_flatsharp_dd0_hdg2 = CuSparseMatrixCSC{Float64}(flatsharp * dual_d0 * hdg_2);
cu_d0 = CuSparseMatrixCSC{Float64}(d0);

cu_codif_1 = CuSparseMatrixCSC{Float64}(codif_1);

lap_first_term = CuSparseMatrixCSC{Float64}(inv_hdg_1 * dual_d0 * hdg_2 * d1);
lap_second_term = CuSparseMatrixCSC{Float64}(d0 * inv_hdg_0 * dual_d1 * hdg_1);
momentum_one_laplacian = lap_first_term + lap_second_term;

function momentum_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, debugat, debug)
  ρ = ρₕ .+ ρ′
  u = wdg_10(U, 1 ./ ρ)

  flow_creation = -wdg_10(U, cu_codif_1 * u) # U ∧ δu
  momentum_advection = -momentum_flatsharp_dd0_hdg2 * wdg_11(u, momentum_flatsharp_hdg_1 * U) - # L(u, U)
                       momentum_flatsharp_hdg_1 * wdg_10(u, momentum_interp_hdg2_d1 * U)
  energy = 0.5 * wdg_10(momentum_flatsharp_dd0_hdg2 * wdg_11(u, momentum_flatsharp_hdg_1 * u), ρ) # 1/2 * ρ * d||u||^2
  pressure_diff = -cu_d0 * perturbed_pressure(Thetaₕ, Theta′) # dP′
  momentum_diff = μ * momentum_one_laplacian * u # μΔu
  body_forces = wdg_10(g, ρ′) # ρ′g

  if debug && step % debugat == 0
    println("##### MOMENTUM DEBUG AT STEP $(step) #####")
    println("Flow Creation: $(CUDA.extrema(flow_creation))")
    println("Momentum Advection: $(CUDA.extrema(momentum_advection))")
    println("Energy: $(CUDA.extrema(energy))")
    println("Pressure Difference: $(CUDA.extrema(pressure_diff))")
    println("Momentum Diffusion: $(CUDA.extrema(momentum_diff))")
    println("Body Forces: $(CUDA.extrema(body_forces))")
    println("##### MOMENTUM DEBUG END #####")
  end

  res = CUDA.zeros(Float64, ne(sd))
  res .= flow_creation .+ momentum_advection .+ energy .+ pressure_diff .+ momentum_diff .+ body_forces

  return res
end

# TODO: Below might have better interpolation but dual_d0 has problems at the boundary
# ptemp_flatsharp_hdg1_d0 = CuSparseMatrixCSC{Float64}(inv_hdg_1 * dual_d0 * dual_interp)
ptemp_flatsharp_hdg1_d0 = CuSparseMatrixCSC{Float64}(flatsharp * hdg_1 * d0);
ptemp_interp_hdg2 = CuSparseMatrixCSC{Float64}(interp * hdg_2);
ptemp_zero_laplacian = CuSparseMatrixCSC{Float64}(inv_hdg_0 * dual_d1 * hdg_1 * d0);

function potential_temperature_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, debugat, debug)
  ρ = ρₕ .+ ρ′
  Theta = Thetaₕ .+ Theta′

  u = wdg_10(U, 1 ./ ρ)
  theta = Theta ./ ρ

  temperature_creation = -Theta .* (cu_codif_1 * u) # Theta ∧ δu
  temperature_advection = -ptemp_interp_hdg2 * wdg_11(u, ptemp_flatsharp_hdg1_d0 * Theta) # L(u, Theta)
  temperature_diffusion = κ * ptemp_zero_laplacian * theta # κΔtheta

  if debug && step % debugat == 0
    println("##### POTENTIAL TEMPERATURE DEBUG AT STEP $(step) #####")
    println("Temperature Creation: $(CUDA.extrema(temperature_creation))")
    println("Temperature Advection: $(CUDA.extrema(temperature_advection))")
    println("Temperature Diffusion: $(CUDA.extrema(temperature_diffusion))")
    println("##### POTENTIAL TEMPERATURE DEBUG END #####")
  end

  res = CUDA.zeros(Float64, nv(sd))
  res .= temperature_creation .+ temperature_advection .+ temperature_diffusion

  return res
end

function run_compressible_ns(Thetaₕ, Theta′₀, U₀, ρₕ, ρ′₀, tₑ, Δt; saveat=500, debug=false)

  Theta′ = deepcopy(Theta′₀)
  U = deepcopy(U₀)
  ρ′ = deepcopy(ρ′₀)

  m₀ = sum(ρₕ .+ ρ′₀)
  E₀ = sum(Thetaₕ .+ Theta′₀)

  Theta′_half = CUDA.zeros(Float64, nv(sd))

  U_half = CUDA.zeros(Float64, ne(sd))

  ρ′_half = CUDA.zeros(Float64, nv(sd))
  ρ′_full = CUDA.zeros(Float64, nv(sd))

  steps = ceil(Int64, tₑ / Δt)

  Us = [Array(U₀)]
  Thetas = [Array(Theta′₀)]
  rhos = [Array(ρ′₀)]

  for step in 1:steps

    apply_periodic!(U, elr_bounds)
    apply_periodic!(ρ′, vlr_bounds)
    apply_periodic!(Theta′, vlr_bounds)

    # apply_periodic!(U, etb_bounds)
    # apply_periodic!(ρ, vtb_bounds)
    # apply_periodic!(Theta, vtb_bounds)

    U_half .= U .+ 0.5 .* Δt * momentum_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, saveat, debug)
    Theta′_half .= Theta′ .+ 0.5 .* Δt * potential_temperature_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, saveat, debug)

    ρ′_full .= rho_smoothing_mat * (ρ′ .- Δt * cu_codif_1 * U_half)
    ρ′_half .= 0.5 .* (ρ′ .+ ρ′_full)

    U .= U .+ Δt * momentum_continuity(Thetaₕ, Theta′_half, U_half, ρₕ, ρ′_half, step, saveat, debug)
    Theta′ .= Theta_smoothing_mat * (Theta′ .+ Δt * potential_temperature_continuity(Thetaₕ, Theta′_half, U_half, ρₕ, ρ′_half, step, saveat, debug))
    ρ′ .= ρ′_full

    if any(CUDA.isnan.(U))
      println("Warning, NAN result in U at step: $(step)")
      break
    elseif any(CUDA.isinf.(U))
      println("Warning, INF result in U at step: $(step)")
      break
    elseif any(CUDA.isnan.(ρ′))
      println("Warning, NAN result in ρ′ at step: $(step)")
      break
    elseif any(CUDA.isinf.(ρ′))
      println("Warning, INF result in ρ′ at step: $(step)")
      break
    elseif any(CUDA.isnan.(Theta′))
      println("Warning, NAN result in Theta′ at step: $(step)")
      break
    elseif any(CUDA.isinf.(Theta′))
      println("Warning, INF result in Theta′ at step: $(step)")
      break
    end

    if step % saveat == 0
      push!(Us, Array(U))
      push!(Thetas, Array(Theta′))
      push!(rhos, Array(ρ′))
      println("Loading simulation results: $((step / steps) * 100)%")
      println("Relative mass is : $((CUDA.sum(ρₕ .+ ρ′) / m₀) * 100)%")
      println("Relative energy is : $((CUDA.sum(Thetaₕ .+ Theta′) / E₀) * 100)%")
      println("-----")

      if step % (5 * saveat) == 0
        fig = Figure()
        ax = CairoMakie.Axis(fig[1, 1]; title="Density Perturbation at Time $(step * Δt)")
        msh = CairoMakie.mesh!(ax, s, color=Array(ρ′), colormap=Reverse(:oslo))
        Colorbar(fig[1, 2], msh)
        display(fig)

        # fig = Figure()
        # ax = CairoMakie.Axis(fig[1, 1]; title="Theta Perturbation at Time $(step * Δt)")
        # msh = CairoMakie.mesh!(ax, s, color=Array(Theta′), colormap=:jet)
        # Colorbar(fig[1, 2], msh)
        # display(fig)

      end

    end
  end

  return Thetas, Us, rhos
end

# For dry air
const Pr = 0.7
Re = Inf
μ = 1 / Re # Momentum diffusivity, 1/Re
κ = μ / Pr # Heat diffusivity

tₑ = 0.1 # 2 # 0.015
Δt = 1e-5 # 2e-5 - 0.04 | 5e-5 - 0.08 | 1e-4 - 0.16 # dx / 360
Thetas, Us, rhos = run_compressible_ns(Thetaₕ, Theta′, U₀, ρₕ, ρ′, tₑ, Δt; saveat=250, debug=false);

filename = "$(simname)_Re=$(Re)_tend=$(tₑ)_dx=$(dx)_rho_$(rho_smooth)_smoothing"
save("$(filename).jld2", Dict("Thetas" => Thetas, "Us" => Us, "rhos" => rhos))

function load_savestate(filename)
  file = load("$(joinpath(pwd(), filename)).jld2")
  return file["Thetas"], file["Us"], file["rhos"]
end

Thetas, Us, rhos = load_savestate(filename);

function reload_initial(Thetas, Us, rhos)
  Theta′ .= CuVector(Thetas[end])
  U₀ .= CuVector(Us[end])
  ρ′ .= CuVector(rhos[end])
end

# reload_initial(Thetas, Us, rhos)

cpu_wdg_10 = dec_wedge_product(Tuple{1,0}, sd);
cpu_wdg_11 = dec_wedge_product(Tuple{1,1}, sd);

cpu_ρₕ = Array(ρₕ)
cpu_Thetaₕ = Array(Thetaₕ)

function view_periodic(x, val = mean(x); tb=false, lr=true)
  y = deepcopy(x)
  if tb
    y[vtb_bounds.boundary_zone] .= val
  end
  if lr
    y[vlr_bounds.boundary_zone] .= val
  end
  return y
end

timestep = length(Us)
u_end = cpu_wdg_10(Us[timestep], 1 ./ (rhos[timestep] .+ cpu_ρₕ))
ω_end = hdg_2 * d1 * u_end

# Velocity Plots
plot_zeroform(s, ω_end; title="Vorticity", colormap=:viridis)
plot_zeroform(s, view_periodic(codif_1 * u_end); title="Velocity Divergence")
plot_zeroform(s, norm.(pp♯ * u_end); title="Velocity Magnitude")
plot_primal_vector_field(s, u_end)

# Momentum Plots
plot_zeroform(s, hdg_2 * d1 * Us[timestep]; title="Density-Coupled Vorticity")
plot_zeroform(s, view_periodic(codif_1 * Us[timestep]); title="Momentum Divergence")
plot_zeroform(s, norm.(pp♯ * Us[timestep]); title="Momentum Magnitude")

# Density Plots
plot_zeroform(s, rhos[timestep] .+ cpu_ρₕ; title="Density", colormap=Reverse(:oslo))
plot_zeroform(s, rhos[timestep]; title="Density Perturbation", colormap=Reverse(:oslo))

# Density-Coupled Potential Temperature Plots

plot_zeroform(s, Thetas[timestep] .+ cpu_Thetaₕ; title="Density-Coupled Potential Temperature")
plot_zeroform(s, Thetas[timestep]; title="Density-Coupled Potential Temperature Perturbation")

# Potential Temperature Plots
theta_end = (Thetas[timestep] .+ cpu_Thetaₕ) ./ (rhos[timestep] .+ cpu_ρₕ)
plot_zeroform(s, theta_end; title="Potential Temperature", colormap=:thermal)
plot_zeroform(s, κ * SparseMatrixCSC(ptemp_zero_laplacian) * theta_end; title="Thermal Diffusion", colormap=:thermal)

# Pressure Plots
plot_zeroform(s, exact_pressure(Thetas[timestep] .+ cpu_Thetaₕ); title="Pressure Field", colormap = Reverse(:acton))
plot_zeroform(s, perturbed_pressure(cpu_Thetaₕ, Thetas[timestep]); title="Pressure Field Perturbation", colormap = Reverse(:acton))

function save_dynamics(save_file_name, endstep, step=1; tb=false, lr=false)
  time = Observable(1)

  U_mag = @lift(norm.(pp♯ * Us[$time]))
  U_vort = @lift(interp * hdg_2 * d1 * cpu_wdg_10(Us[$time], 1 ./ (rhos[$time] .+ cpu_ρₕ)))
  ρ = @lift(view_periodic(rhos[$time], mean(rhos[$time]); tb=tb, lr=lr))
  # theta = @lift(view_periodic(Thetas[$time], mean(Thetas[$time]); tb=tb, lr=lr))
  # P = @lift(perturbed_pressure(cpu_Thetaₕ, view_periodic(Thetas[$time], mean(Thetas[$time]); tb=tb, lr=lr)))
  # U_div = @lift(Array(codif_1 * Us[$time]))
  theta = @lift(view_periodic((Thetas[$time] .+ cpu_Thetaₕ) ./ (rhos[$time] .+ cpu_ρₕ); tb=tb, lr=lr))


  f = Figure()

  ax_U_vort = CairoMakie.Axis(f[1, 1], title="Vorticity")
  msh_U_vort = mesh!(ax_U_vort, s; color=U_vort, colormap=:viridis)
  Colorbar(f[1, 2], msh_U_vort)

  ax_ρ = CairoMakie.Axis(f[2, 1], title="Density Perturbation")
  msh_ρ = mesh!(ax_ρ, s; color=ρ, colormap=Reverse(:oslo))
  Colorbar(f[2, 2], msh_ρ)

  # ax_theta = CairoMakie.Axis(f[1, 3], title="Momentum Divergence")
  # msh_theta = mesh!(ax_theta, s; color=U_div, colormap=:thermal)
  # Colorbar(f[1, 4], msh_theta)

  # ax_theta = CairoMakie.Axis(f[1, 3], title="Density-Coupled Potential Temperature Perturbation")
  # msh_theta = mesh!(ax_theta, s; color=theta, colormap=:thermal)
  # Colorbar(f[1, 4], msh_theta)

  ax_U_div = CairoMakie.Axis(f[1, 3], title="Momentum Magnitude")
  msh_U_div = mesh!(ax_U_div, s; color=U_mag, colormap=:viridis)
  Colorbar(f[1, 4], msh_U_div)

  ax_theta = CairoMakie.Axis(f[2, 3], title="Potential Temperature")
  msh_theta = mesh!(ax_theta, s; color=theta, colormap=:thermal)
  Colorbar(f[2, 4], msh_theta)


  # ax_P = CairoMakie.Axis(f[2, 3], title="Pressure Field Perturbation")
  # msh_P = mesh!(ax_P, s; color=P, colormap=Reverse(:acton))
  # Colorbar(f[2, 4], msh_P)

  timestamps = 1:step:endstep
  CairoMakie.record(f, save_file_name, timestamps; framerate=15) do t
    time[] = t
  end
end

save_dynamics("$(filename).mp4", length(Us), 40; lr=true)



# m₀ = CUDA.sum(ρ₀)
# E₀ = CUDA.sum(Theta₀)

# fig = Figure();
# ax = CairoMakie.Axis(fig[1, 1]; title="Relative Error in Mass", xlabel="Index of Saved Data", ylabel="Relative Error")
# rho_data = (sum.(rhos) .- m₀) ./ m₀
# CairoMakie.plot!(ax, rho_data)
# fig

# fig = Figure();
# ax = CairoMakie.Axis(fig[1, 1]; title="Relative Error in Thermal Energy", xlabel="Index of Saved Data", ylabel="Relative Error")
# Theta_data = (sum.(Thetas) .- E₀) ./ E₀
# CairoMakie.plot!(ax, Theta_data)
# fig
