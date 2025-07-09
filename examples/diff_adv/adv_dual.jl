using CombinatorialSpaces
using CairoMakie
using Distributions
using StaticArrays
using SparseArrays

lx = 150
ly = 60
dx = dy = 1
s = triangulated_grid(lx, ly, dx, dy, Point3d, false)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

dual_d0 = dec_dual_derivative(0, sd);
dual_d1 = dec_dual_derivative(1, sd);

hdg_1 = dec_hodge_star(1, sd, DiagonalHodge())
hdg_2 = dec_hodge_star(2, sd, DiagonalHodge())

inv_hdg_0 = dec_inv_hodge_star(0, sd, DiagonalHodge())
inv_hdg_1 = dec_inv_hodge_star(1, sd, DiagonalHodge())

wdg_10 = dec_wedge_product(Tuple{1,0}, sd)
wdg_11 = dec_wedge_product(Tuple{1,1}, sd)

dd♯ = ♯_mat(sd, LLSDDSharp())

u = eval_constant_primal_form(sd, Point3d(1, 0, 0))

dual_points = sd[triangle_center(sd), :dual_point]

v_dist = MvNormal([lx / 4.0, ly / 2.0], [3, 3])
v♯ = [100 * Point3d(1, 1, 0) * pdf(v_dist, [p[1], p[2]]) for p in dual_points]
v♯ = 0.01 * map(_ -> Point3d(1, 0, 0), dual_points)

x = map(p -> p[1], dual_points)
y = map(p -> p[2], dual_points)

vel_u = map(p -> p[1], v♯)
vel_v = map(p -> p[2], v♯)

fig = Figure(size=(1000, 1000));
ax = CairoMakie.Axis(fig[1, 1])
CairoMakie.arrows!(ax, x, y, vel_u, vel_v)
# CairoMakie.wireframe!(ax, sd)
fig

function dd♭_mat(sd)
  mat_type = SMatrix{1,length(eltype(s[:point])),eltype(eltype(s[:point])),length(eltype(s[:point]))}
  ♭_mat = spzeros(mat_type, ne(s), ntriangles(s))
  for e in edges(s)
    des = elementary_duals(1, sd, e)
    for de in des
      dv = first(dual_edge_vertices(sd, de))
      tgt = sd[de, :D_∂v0]
      src = sd[de, :D_∂v1]
      de_sign = sd[de, :D_edge_orientation] ? 1 : -1
      de_vec = (dual_point(sd, tgt) - dual_point(sd, src)) * de_sign

      ♭_mat[e, dv-nv(sd)-ne(sd)] = mat_type(de_vec)
    end
  end
  ♭_mat
end

dd♭ = dd♭_mat(sd)

v₀ = only.(dd♭ * v♯)

test = dd♯ * v₀

test_u = map(p -> p[1], test)
test_v = map(p -> p[2], test)

fig = Figure(size=(1000, 1000));
ax = CairoMakie.Axis(fig[1, 1])
CairoMakie.arrows!(ax, x, y, test_u, test_v)
fig

form_zero_interp = inv_hdg_0 * dual_d1
form_two_interp = dual_d0 * hdg_2

Δt = 1e-2
tₑ = 80
v = deepcopy(v₀)

for t in 1:ceil(Int64, tₑ / Δt)
  v .= v .- Δt .* (form_two_interp * wdg_11(u, inv_hdg_1 * v) + hdg_1 * wdg_10(u, form_zero_interp * v))
end

res = dd♯ * v
res_u = map(p -> p[1], res)
res_v = map(p -> p[2], res)

fig = Figure(size=(1000, 1000));
ax = CairoMakie.Axis(fig[1, 1])
CairoMakie.arrows!(ax, x, y, res_u, res_v)
fig
