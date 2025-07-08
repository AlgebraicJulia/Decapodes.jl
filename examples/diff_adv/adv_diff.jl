using CairoMakie
using CombinatorialSpaces
using ComponentArrays
using Distributions
using LinearAlgebra
using MLStyle
using SparseArrays
using StaticArrays

lx = ly = 1
dx = dy = 0.01
s = triangulated_grid(lx, ly, dx, dy, Point3d, false)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
CairoMakie.wireframe!(ax, s)
fig

d0 = dec_differential(0, sd)

dual_d0 = dec_dual_derivative(0, sd);
dual_d1 = dec_dual_derivative(1, sd);

# dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));

hdg_1 = dec_hodge_star(1, sd, DiagonalHodge())
hdg_2 = dec_hodge_star(2, sd, DiagonalHodge())

inv_hdg_0 = dec_inv_hodge_star(0, sd, DiagonalHodge())
inv_hdg_2 = dec_inv_hodge_star(2, sd, DiagonalHodge())

wdg_11 = dec_wedge_product(Tuple{1,1}, sd)

interp = SparseMatrixCSC{Float64}(inv_hdg_0) * p2_d2_interpolation(sd) * SparseMatrixCSC{Float64}(inv_hdg_2)
flatsharp = ♭♯_mat(sd)

ptemp_flatsharp_hdg1_d0 = flatsharp * hdg_1 * d0
ptemp_interp_hdg2 = interp * hdg_2
ptemp_zero_laplacian = inv_hdg_0 * dual_d1 * hdg_1 * d0

# function select_left_vertical_bar(dx, lx, s, idx)
#     nx = length(0:dx:lx)
#     shift = idx < 0 ? nx + idx : idx
#     collect(1:nx:nv(s)) .+ shift
# end

# function select_right_vertical_bar(dx, lx, s, idx)
#     nx = length(0:dx:lx)
#     shift = idx < 0 ? nx + idx - 1 : idx
#     lw = collect(1:nx:nv(s)) .+ shift
#     for (i, _) in enumerate(lw)
#         if isodd(i)
#             lw[i] += 1
#         end
#     end
#     return lw
# end

# function select_top_horizontal_bar(dx, lx, s, idx)
#     nx = length(0:dx:lx)
#     shift = idx < 0 ? nx + idx - 1 : idx
#     lw = collect(1:2nx:nv(s)) .+ shift
#     return vcat(lw, lw .+ 1)
# end

# function select_bottom_horizontal_bar(dx, lx, s, idx)
#     nx = length(0:dx:lx)
#     shift = idx < 0 ? nx + idx - 1 : idx
#     lw = collect(nx+1:2nx:nv(s)) .+ shift
#     return vcat(lw, lw .+ 1)
# end

# function find_edge(s, p0, p1)
#     for e in edges(s)
#         if (p0 in s[e, :∂v0] || p0 in s[e, :∂v1]) &&
#            (p1 in s[e, :∂v0] || p1 in s[e, :∂v1])
#             return e
#         end
#     end
#     return 0
# end

# function get_vertical_edges(s, arr)
#     es = Int64[]
#     for (a1, a2) in zip(arr, arr[2:end])
#         push!(es, find_edge(s, a1, a2))
#     end
#     return es
# end

# function get_horizontal_edges(s, arr)
#     es = Int64[]
#     for (a1, a2) in zip(arr, arr[length(arr)÷2+1:end])
#         push!(es, find_edge(s, a1, a2))
#     end
#     return es
# end

# left_boundary_zone_one_form = vcat(
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, 1)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, 1)),
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, 0)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, 0)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, 0)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, 1)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, 0)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, 1)))

# copy_to_left_zone_one_form = vcat(
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, -3)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, -4)),
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, -4)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, -5)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, -4)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, -3)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, -4)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, -3)))

# copy_to_right_zone_one_form = vcat(
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, 2)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, 2)),
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, 3)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, 3)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, 2)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, 3)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, 2)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, 3)))

# right_boundary_zone_one_form = vcat(
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, -2)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, -2)),
#     get_vertical_edges(s, select_right_vertical_bar(dx, lx, s, -1)),
#     get_vertical_edges(s, select_left_vertical_bar(dx, lx, s, -1)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, -2)),
#     get_horizontal_edges(s, select_top_horizontal_bar(dx, lx, s, -1)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, -2)),
#     get_horizontal_edges(s, select_bottom_horizontal_bar(dx, lx, s, -1)))

# left_boundary_zone_zero_form = vcat(
#     select_right_vertical_bar(dx, lx, s, 1),
#     select_left_vertical_bar(dx, lx, s, 1),
#     select_right_vertical_bar(dx, lx, s, 0),
#     select_left_vertical_bar(dx, lx, s, 0),
#     select_top_horizontal_bar(dx, lx, s, 0),
#     select_top_horizontal_bar(dx, lx, s, 1),
#     select_bottom_horizontal_bar(dx, lx, s, 0),
#     select_bottom_horizontal_bar(dx, lx, s, 1))

# copy_to_left_zone_zero_form = vcat(
#     select_right_vertical_bar(dx, lx, s, -3),
#     select_left_vertical_bar(dx, lx, s, -4),
#     select_right_vertical_bar(dx, lx, s, -4),
#     select_left_vertical_bar(dx, lx, s, -5),
#     select_top_horizontal_bar(dx, lx, s, -4),
#     select_top_horizontal_bar(dx, lx, s, -3),
#     select_bottom_horizontal_bar(dx, lx, s, -4),
#     select_bottom_horizontal_bar(dx, lx, s, -3))

# copy_to_right_zone_zero_form = vcat(
#     select_right_vertical_bar(dx, lx, s, 2),
#     select_left_vertical_bar(dx, lx, s, 2),
#     select_right_vertical_bar(dx, lx, s, 3),
#     select_left_vertical_bar(dx, lx, s, 3),
#     select_top_horizontal_bar(dx, lx, s, 2),
#     select_top_horizontal_bar(dx, lx, s, 3),
#     select_bottom_horizontal_bar(dx, lx, s, 2),
#     select_bottom_horizontal_bar(dx, lx, s, 3))

# right_boundary_zone_zero_form = vcat(
#     select_right_vertical_bar(dx, lx, s, -2),
#     select_left_vertical_bar(dx, lx, s, -2),
#     select_right_vertical_bar(dx, lx, s, -1),
#     select_left_vertical_bar(dx, lx, s, -1),
#     select_top_horizontal_bar(dx, lx, s, -2),
#     select_top_horizontal_bar(dx, lx, s, -1),
#     select_bottom_horizontal_bar(dx, lx, s, -2),
#     select_bottom_horizontal_bar(dx, lx, s, -1))

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

# e_bounds = EdgeMapping(sd, bounds)
begin
    # v = 2
    # p1 = edge_vertices(sd, e_bounds.boundary_zone)[v]
    # p2 = edge_vertices(sd, e_bounds.copy_to_boundary_zone)[v]

    # p1 = vtb_bounds.boundary_zone
    # p2 = vtb_bounds.copy_to_boundary_zone

    p3 = vlr_bounds.boundary_zone
    p4 = vlr_bounds.copy_to_boundary_zone

    test = zeros(nv(sd))

    # test[p1] .= 1
    # test[p2] .= 2

    test[p3] .= 3
    test[p4] .= 4

    fig = Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    msh = CairoMakie.mesh!(ax, s; color=test)
    CairoMakie.Colorbar(fig[1, 2], msh)
    # CairoMakie.wireframe!(ax, s)
    display(fig)
end

tb_vertex_periodic_boundary = VertexMapping(sd, tb_bounds)
lr_vertex_periodic_boundary = VertexMapping(sd, lr_bounds)
tb_edge_periodic_boundary = EdgeMapping(sd, tb_bounds)
lr_edge_periodic_boundary = EdgeMapping(sd, lr_bounds)

function potential_temperature_continuity(theta, u)
    return -ptemp_interp_hdg2 * wdg_11(u, ptemp_flatsharp_hdg1_d0 * theta) + # L(u, theta)
           κ * ptemp_zero_laplacian * theta # κΔtheta 

end

function run_heat(theta₀, u, tₑ, Δt; saveat=500)

    theta = deepcopy(theta₀)

    steps = ceil(Int64, tₑ / Δt)

    thetas = [deepcopy(theta₀)]

    for step in 1:steps
        apply_periodic!(theta, lr_vertex_periodic_boundary)
        # apply_periodic!(theta, tb_vertex_periodic_boundary)

        theta .= theta .+ Δt * potential_temperature_continuity(theta, u)

        if step % saveat == 0
            push!(thetas, deepcopy(theta))
            println("Loading simulation results: $((step / steps) * 100)%")
            println("-----")
        end
    end

    return thetas
end

κ = 1e-3
tₑ = 5
Δt = 1e-3

u = eval_constant_primal_form(sd, Point3d(5e-1, 0, 0))
# u = eval_constant_primal_form(sd, Point3d(0, -5e-1, 0))
# u = eval_constant_primal_form(sd, Point3d(5e-1, 5e-1, 0))

theta_dist = MvNormal([lx / 2, ly / 2], [0.1, 0.1])
theta_dist_1 = MvNormal([lx / 4, ly * 3 / 4], [0.1, 0.1])
theta_dist_2 = MvNormal([lx * 3 / 4, ly / 4], [0.1, 0.1])

theta = rand(Float64, nv(sd)) .+ [pdf(theta_dist_1, [p[1], p[2]]) for p in sd[:point]] .+ [pdf(theta_dist_2, [p[1], p[2]]) for p in sd[:point]] .+ 300
theta = rand(Float64, nv(sd)) .+ [pdf(theta_dist, [p[1], p[2]]) for p in sd[:point]] .+ 300

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1]; title="Initial Potential Temperature")
msh = CairoMakie.mesh!(ax, s, color=theta, colormap=:jet)
Colorbar(fig[1, 2], msh)
display(fig)

thetas = run_heat(theta, u, tₑ, Δt; saveat=100)

for t in thetas
    fig = Figure()
    ax = CairoMakie.Axis(fig[1, 1]; title="Potential Temperature")
    msh = CairoMakie.mesh!(ax, s, color=t, colormap=:jet)
    Colorbar(fig[1, 2], msh)
    display(fig)
end