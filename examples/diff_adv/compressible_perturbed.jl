using CairoMakie
using CombinatorialSpaces
using ComponentArrays
using Distributions
using LinearAlgebra
using MLStyle
using SparseArrays
using StaticArrays

using CUDA
using CUDA.CUSPARSE

lx = ly = 6
dx = dy = 0.04
s = triangulated_grid(lx, ly, dx, dy, Point3d, false)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
wireframe!(ax, s)
fig

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
flatsharp = ♭♯_mat(sd)

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
topb = collect_top_boundary_triangles(dx, lx, dy, ly, depth)
botb = collect_bottom_boundary_triangles(dx, lx, dy, ly, depth)

leftb = collect_left_boundary_triangles(dx, lx, dy, ly, depth)
rightb = collect_right_boundary_triangles(dx, lx, dy, ly, depth)

import Base.vcat
vcat(b1::AbstractBoundaryMapping, b2::AbstractBoundaryMapping) = TriangleMapping(vcat(b1.boundary_zone, b2.boundary_zone),
  vcat(b1.copy_to_boundary_zone, b2.copy_to_boundary_zone))

function apply_periodic!(x, m::AbstractBoundaryMapping)
  x[m.boundary_zone] .= x[m.copy_to_boundary_zone]
  return x
end

tb_bounds = vcat(topb, botb)
lr_bounds = vcat(leftb, rightb)
vtb_bounds = VertexMapping(sd, tb_bounds)
vlr_bounds = VertexMapping(sd, lr_bounds)
etb_bounds = EdgeMapping(sd, tb_bounds)
elr_bounds = EdgeMapping(sd, lr_bounds)

##############################
##### INITIAL CONDITIONS #####
##############################

P₀ = 1e5 # Reference pressure for potential temperature
Cₚ = 1006 # Specific heat
R = 287 # Specific gas constant
ρᵣ = 1 # Reference density
θᵣ = 300 # Reference potential temperature, temperature at P₀
Pᵣ = ρᵣ * R * θᵣ
gₐ = -20_000 # -9.81 # Acceleration due to gravity
g = CuVector{Float64}(eval_constant_primal_form(sd, Point3d(0, gₐ, 0)))

# TODO: This assumes theta is constant in z
function hydrostatic_pressure(theta, h)
  return (1 / (Cₚ * theta) * (P₀)^(R / Cₚ) * gₐ * h + Pᵣ^(R / Cₚ))^(Cₚ / R)
end

function hydrostatic_density(theta, h)
  Pₕ = hydrostatic_pressure(theta, h)

  return (Pₕ / (R * theta)) * (P₀ / Pₕ)^(R / Cₚ)
end

hydrostatic_pressure(300, 0) ≈ Pᵣ
hydrostatic_density(300, 0)

# # CONSTANT
# simname = "Perturbed_Constant_Test"
# ρ₀ = CUDA.ones(Float64, nv(sd))
# U₀ = CUDA.zeros(Float64, ne(sd))
# theta₀ = 300 * CUDA.ones(Float64, nv(sd))
# Theta₀ = ρ₀ .* theta₀

# HYDROSTATIC
simname = "Perturbed_Hydrostatic"
thetaₕ = 300 * CUDA.ones(Float64, nv(sd))
ρₕ = CuVector{Float64}(map(p -> hydrostatic_density(300, p[2]), sd[:point]))
Thetaₕ = ρₕ .* thetaₕ
U₀ = CUDA.zeros(Float64, ne(sd))
ρ′ = CUDA.zeros(Float64, nv(sd))
Theta′ = CUDA.zeros(Float64, nv(sd))

# # WARM BUBBLE 
# simname = "Perturbed_Warm_Bubble"
# ρ₀ = CUDA.ones(Float64, nv(sd))
# U₀ = CUDA.zeros(Float64, ne(sd))
# theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
# theta₀ = 300 .+ 0.5 * CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])
# Theta₀ = ρ₀ .* theta₀

# WARM BUBBLE IN DENSITY GRADIENT
simname = "Perturbed_Warm_Bubble_Gradient"
thetaₕ = 300
ρ₀ = CuVector{Float64}(map(p -> hydrostatic_density(thetaₕ, p[2]), sd[:point]))
U₀ = CUDA.zeros(Float64, ne(sd))
theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
theta₀ = thetaₕ * CUDA.ones(Float64, nv(sd)) .+ 5 * CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])
Theta₀ = ρ₀ .* theta₀

# # COLD BUBBLE
# simname = "Perturbed_Cold_Bubble"
# ρ₀ = CUDA.ones(Float64, nv(sd))
# U₀ = CUDA.zeros(Float64, ne(sd))
# theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
# theta₀ = 300 .- 0.05 * CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])
# Theta₀ = ρ₀ .* theta₀

# # COLD BUBBLE WITH WIND
# simname = "Perturbed_Cold_Bubble_Wind"
# ρ₀ = CUDA.ones(Float64, nv(sd))
# U₀ = CuVector{Float64}(eval_constant_primal_form(sd, Point3d(100, 0, 0))) # Techincally velocity but we have unit mass
# theta_dist = MvNormal([lx / 2, ly / 2], 0.1)
# theta₀ = 300 .- 0.05 * CuArray{Float64}([pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]])
# Theta₀ = ρ₀ .* theta₀

# LIGHT BUBBLE IN DENSITY GRADIENT
simname = "Perturbed_Light_Bubble_Gradient"
thetaₕ = 300 * CUDA.ones(Float64, nv(sd))
ρₕ = CuVector{Float64}(map(p -> hydrostatic_density(300, p[2]), sd[:point]))
Thetaₕ = ρₕ .* thetaₕ

U₀ = CUDA.zeros(Float64, ne(sd))
rho_dist = MvNormal([lx / 2, ly / 4], 0.1)
ρ′ = -0.01 * CuArray{Float64}([pdf(rho_dist, [p[1], p[2]]) for p in sd[:point]])
Theta′ = CUDA.zeros(Float64, nv(sd))

# # RAYLEIGH-TAYLOR INSTABILITY
# simname = "Perturbed_RTI"
# ρ₀ = CUDA.ones(Float64, nv(sd))
# U₀ = CUDA.zeros(Float64, ne(sd))
# theta₀ = CUDA.zeros(Float64, nv(sd))

# upper_points = findall(p -> p[2] >= lx * 0.75, sd[:point])
# lower_points = findall(p -> p[2] < lx * 0.75, sd[:point])

# theta₀[upper_points] .= 302
# theta₀[lower_points] .= 300
# theta₀ .+= 1e-6 * ((2 * CUDA.rand(Float64, nv(sd))) .- 1)
# Theta₀ = ρ₀ .* theta₀

##################################
##### END INITIAL CONDITIONS #####
##################################

boundary_edges = findall(x -> x != 0, dual_d0 * ones(ntriangles(sd)))
boundary_vertices = unique(vcat(sd[boundary_edges, :∂v0], sd[boundary_edges, :∂v1]))

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

plot_zeroform(s, Array(ρₕ))
plot_zeroform(s, Array(Thetaₕ))
plot_zeroform(s, Array(exact_pressure(Thetaₕ)))

plot_zeroform(s, Array(ρ′))
plot_zeroform(s, Array(Theta′))
plot_zeroform(s, Array(perturbed_pressure(Thetaₕ, Theta′)))

function smoothing(sd, c_smooth)
  mat = spzeros(nv(sd), nv(sd))
  for e in edges(sd)
    v1 = sd[e, :∂v0]
    v2 = sd[e, :∂v1]
    w = 1 / sd[e, :length]
    mat[v1, v2] = w
    mat[v2, v1] = w
  end

  c = c_smooth ./ 2

  for v in vertices(sd)
    row = mat[v, :]
    tot_w = sum(row)
    for i in row.nzind
      mat[v, i] = c .* row[i] ./ tot_w
    end

    mat[v, v] = (1 - c)
  end
  return mat
end

c_smooth = 0.2
forward_smooth = smoothing(sd, c_smooth)
backward_smooth = smoothing(sd, -c_smooth)

smoothing_mat = CuSparseMatrixCSC{Float64}(backward_smooth * forward_smooth)

momentum_flatsharp_hdg_1 = CuSparseMatrixCSC{Float64}(flatsharp * hdg_1)
momentum_interp_hdg2_d1 = CuSparseMatrixCSC{Float64}(interp * abs.(hdg_2) * d1)
momentum_flatsharp_dd0_hdg2 = CuSparseMatrixCSC{Float64}(flatsharp * dual_d0 * hdg_2)
cu_d0 = CuSparseMatrixCSC{Float64}(d0)

codif_1 = CuSparseMatrixCSC{Float64}(SparseMatrixCSC(inv_hdg_0) * dual_d1 * SparseMatrixCSC(hdg_1))

lap_first_term = CuSparseMatrixCSC{Float64}(inv_hdg_1 * dual_d0 * abs.(hdg_2) * d1)
lap_second_term = CuSparseMatrixCSC{Float64}(d0 * inv_hdg_0 * dual_d1 * hdg_1)
momentum_one_laplacian = lap_first_term + lap_second_term

function momentum_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, debugat, debug)
  ρ = ρₕ .+ ρ′
  u = wdg_10(U, 1 ./ ρ)

  flow_creation = -wdg_10(U, codif_1 * u) # U ∧ δu
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

  return flow_creation + momentum_advection + energy + pressure_diff + momentum_diff + body_forces

end

# TODO: Below might have better interpolation but dual_d0 has problems at the boundary
# ptemp_flatsharp_hdg1_d0 = CuSparseMatrixCSC{Float64}(inv_hdg_1 * dual_d0 * dual_interp)
ptemp_flatsharp_hdg1_d0 = CuSparseMatrixCSC{Float64}(flatsharp * hdg_1 * d0)
ptemp_interp_hdg2 = CuSparseMatrixCSC{Float64}(interp * hdg_2)
ptemp_zero_laplacian = CuSparseMatrixCSC{Float64}(inv_hdg_0 * dual_d1 * hdg_1 * d0)

function potential_temperature_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, debugat, debug)
  ρ = ρₕ .+ ρ′
  Theta = Thetaₕ .+ Theta′

  u = wdg_10(U, 1 ./ ρ)
  theta = Theta ./ ρ

  temperature_creation = -Theta .* (codif_1 * u) # Theta ∧ δu
  temperature_advection = -ptemp_interp_hdg2 * wdg_11(u, ptemp_flatsharp_hdg1_d0 * Theta) # L(u, Theta)
  temperature_diffusion = κ * ptemp_zero_laplacian * theta # κΔtheta

  if debug && step % debugat == 0
    println("##### POTENTIAL TEMPERATURE DEBUG AT STEP $(step) #####")
    println("Temperature Creation: $(CUDA.extrema(temperature_creation))")
    println("Temperature Advection: $(CUDA.extrema(temperature_advection))")
    println("Temperature Diffusion: $(CUDA.extrema(temperature_diffusion))")
    println("##### POTENTIAL TEMPERATURE DEBUG END #####")
  end

  return temperature_creation + temperature_advection + temperature_diffusion

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

  Us = [deepcopy(U₀)]
  Thetas = [deepcopy(Theta′₀)]
  rhos = [deepcopy(ρ′₀)]

  for step in 1:steps

    apply_periodic!(U, elr_bounds)
    apply_periodic!(ρ′, vlr_bounds)
    apply_periodic!(Theta′, vlr_bounds)

    # apply_periodic!(U, etb_bounds)
    # apply_periodic!(ρ, vtb_bounds)
    # apply_periodic!(Theta, vtb_bounds)

    U_half .= U .+ 0.5 .* Δt * momentum_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, saveat, debug)
    Theta′_half .= Theta′ .+ 0.5 .* Δt * potential_temperature_continuity(Thetaₕ, Theta′, U, ρₕ, ρ′, step, saveat, debug)

    ρ′_full .= smoothing_mat * (ρ′ .- Δt * codif_1 * U_half)
    ρ′_half .= 0.5 .* (ρ′ .+ ρ′_full)

    U .= U .+ Δt * momentum_continuity(Thetaₕ, Theta′_half, U_half, ρₕ, ρ′_half, step, saveat, debug)
    Theta′ .= Theta′ .+ Δt * potential_temperature_continuity(Thetaₕ, Theta′_half, U_half, ρₕ, ρ′_half, step, saveat, debug)
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
      push!(Us, deepcopy(U))
      push!(Thetas, deepcopy(Theta′))
      push!(rhos, deepcopy(ρ′))
      println("Loading simulation results: $((step / steps) * 100)%")
      println("Relative mass is : $((CUDA.sum(ρₕ .+ ρ′) / m₀) * 100)%")
      println("Relative energy is : $((CUDA.sum(Thetaₕ .+ Theta′) / E₀) * 100)%")
      println("-----")
    end
  end

  return Thetas, Us, rhos
end

# For dry air
Pr = 0.7
μ = 1e-3 # Momentum diffusivity, 1/Re
κ = μ / Pr # Heat diffusivity

tₑ = 0.5
Δt = 1e-5 # dx / 360
Thetas, Us, rhos = run_compressible_ns(Thetaₕ, Theta′, U₀, ρₕ, ρ′, tₑ, Δt; saveat=100, debug=false);

timestep = length(Us)
u_end = wdg_10(Us[timestep], 1 ./ (rhos[timestep] .+ ρₕ))
ω_end = abs.(hdg_2) * d1 * Array(u_end)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Vorticity")
msh = CairoMakie.mesh!(ax, s, color=interp * ω_end, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Velocity Divergence")
msh = CairoMakie.mesh!(ax, s, color=Array(codif_1 * u_end), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Velocity Magnitude")
msh = CairoMakie.mesh!(ax, s, color=norm.(pp♯ * Array(u_end)), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

# Velocity components
plot_primal_vector_field(s, Array(u_end))

ρω_end = abs.(hdg_2) * d1 * Array(Us[timestep])

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Density-Coupled Vorticity")
msh = CairoMakie.mesh!(ax, s, color=interp * ρω_end, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Momentum Divergence")
msh = CairoMakie.mesh!(ax, s, color=Array(codif_1 * Us[timestep]), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Momentum Magnitude")
msh = CairoMakie.mesh!(ax, s, color=norm.(pp♯ * Array(Us[timestep])), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Density")
msh = CairoMakie.mesh!(ax, s, color=Array(rhos[timestep] .+ ρₕ), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Density Perturbation")
msh = CairoMakie.mesh!(ax, s, color=Array(rhos[timestep]), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Density-Coupled Potential Temperature")
msh = CairoMakie.mesh!(ax, s, color=Array(Thetas[timestep] .+ Thetaₕ), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Density-Coupled Potential Temperature Perturbation")
msh = CairoMakie.mesh!(ax, s, color=Array(Thetas[timestep]), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

theta_end = (Thetas[timestep] .+ Thetaₕ) ./ (rhos[timestep] .+ ρₕ)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Potential Temperature")
msh = CairoMakie.mesh!(ax, s, color=Array(theta_end), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Pressure Field")
msh = CairoMakie.mesh!(ax, s, color=Array(exact_pressure(Thetas[timestep] .+ Thetaₕ)), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Pressure Field Perturbation")
msh = CairoMakie.mesh!(ax, s, color=Array(perturbed_pressure(Thetaₕ, Thetas[timestep])), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Thermal Diffusion")
msh = CairoMakie.mesh!(ax, s, color=Array(ptemp_zero_laplacian * theta_end), colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

function view_periodic(x, val; tb=true, lr=true)
  y = deepcopy(x)
  if tb
    y[vtb_bounds.boundary_zone] .= val
  end
  if lr
    y[vlr_bounds.boundary_zone] .= val
  end
  return y
end

function save_dynamics(save_file_name, endstep, step=1; tb=false, lr=false)
  time = Observable(1)

  U_mag = @lift(norm.(pp♯ * Array(Us[$time])))
  ρ = @lift(Array(view_periodic(rhos[$time], mean(rhos[$time]); tb=tb, lr=lr)))
  theta = @lift(Array(view_periodic(Thetas[$time], mean(Thetas[$time]); tb=tb, lr=lr)))
  P = @lift(Array(exact_pressure(view_periodic(Thetas[$time], mean(Thetas[$time]); tb=tb, lr=lr))))
  # U_div = @lift(Array(codif_1 * Us[$time]))

  f = Figure()

  ax_U_div = CairoMakie.Axis(f[1, 1], title="Momentum Magnitude")
  msh_U_div = mesh!(ax_U_div, s; color=U_mag, colormap=:viridis)
  Colorbar(f[1, 2], msh_U_div)

  ax_ρ = CairoMakie.Axis(f[2, 1], title="Density Perturbation")
  msh_ρ = mesh!(ax_ρ, s; color=ρ, colormap=Reverse(:oslo))
  Colorbar(f[2, 2], msh_ρ)

  # ax_theta = CairoMakie.Axis(f[1, 3], title="Momentum Divergence")
  # msh_theta = mesh!(ax_theta, s; color=U_div, colormap=:thermal)
  # Colorbar(f[1, 4], msh_theta)

  ax_theta = CairoMakie.Axis(f[1, 3], title="Density-Coupled Potential Temperature Perturbation")
  msh_theta = mesh!(ax_theta, s; color=theta, colormap=:thermal)
  Colorbar(f[1, 4], msh_theta)

  ax_P = CairoMakie.Axis(f[2, 3], title="Pressure Field Perturbation")
  msh_P = mesh!(ax_P, s; color=P, colormap=Reverse(:acton))
  Colorbar(f[2, 4], msh_P)

  timestamps = 1:step:endstep
  CairoMakie.record(f, save_file_name, timestamps; framerate=15) do t
    time[] = t
  end
end

save_dynamics("$(simname)_mu=$(μ).mp4", length(Us), 10; lr=true)

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