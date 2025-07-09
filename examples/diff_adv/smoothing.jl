using CombinatorialSpaces
using CairoMakie
using SparseArrays

lx = ly = 10
dx = dy = 0.05

s = triangulated_grid(lx, ly, dx, dy, Point3d)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s);
subdivide_duals!(sd, Circumcenter());

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
wireframe!(ax, s)
fig

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

c_smooth = 1
forward_smooth = smoothing(sd, c_smooth)
backward_smooth = smoothing(sd, -c_smooth)

test = map(p -> sin(p[1]) + 0.1 * sin(5 * p[1]) + 0.01 * sin(10 * p[1]), sd[:point])

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s; color=test, colorrange=extrema(test))
Colorbar(fig[1, 2], msh)
fig

res = backward_smooth * forward_smooth * test

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s; color=res, colorrange=extrema(test))
Colorbar(fig[1, 2], msh)
fig

fig = Figure();
ax = CairoMakie.Axis(fig[1, 1])
msh = CairoMakie.mesh!(ax, s; color=res .- test)
Colorbar(fig[1, 2], msh)
fig
