using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using GLMakie
using Distributions

lx = 100
ly = lz = 50

s = parallelepiped(lx = lx, ly = ly, lz = lz; tetcmd = "vpVq2a2.5")
sd = EmbeddedDeltaDualComplex3D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

# GLMakie.wireframe(s)

wdg_01 = dec_wedge_product(Tuple{0,1}, sd)
star₁_mat = hodge_star(1,sd,DiagonalHodge())
dual_d₂_mat = dual_derivative(2,sd)
inv_star₃_mat = inv_hodge_star(0,sd,DiagonalHodge())

codif = inv_star₃_mat * dual_d₂_mat * star₁_mat

T_dist = MvNormal([30, ly/2, lz/2], [10,10,10])
T = [pdf(T_dist, [p[1], p[2], p[3]]) for p in sd[:point]] / 1000
T₀ = T
GLMakie.mesh(s; color = T, transparency=true, shading=NoShading)

c = 1

#TODO: Work on this, same as 2D basically
u = eval_constant_primal_form(sd, Point3d(1,0,0))

u₀ = ComponentArray(T=T, u=u)

function adv(du, u, p, t)
    du.T .-= p.c * inv_star₃_mat * dual_d₂_mat * star₁_mat * wdg_01(u.T, u.u)
end

tₑ = 10
prob = ODEProblem(adv, u₀, (0, tₑ), (c=c,))
soln = solve(prob, Tsit5(), saveat = 0.1)

GLMakie.mesh(s; color = soln.u[end].T, transparency=true, shading=NoShading, colorrange = extrema(T₀))

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)

  f = Figure()

  ls_T = GLMakie.LScene(f[1, 1])
  msh_T = GLMakie.mesh!(ls_T, s; color=T, colorrange = extrema(T₀), transparency=true, shading=NoShading)
  Colorbar(f[1,2], msh_T)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("3d_adv.mp4", 30)