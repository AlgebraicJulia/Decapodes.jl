using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using GLMakie
using Distributions

lx = ly = lz = 10

s = parallelepiped(lx = lx, ly = ly, lz = lz; tetcmd = "vpq1.5a0.05")
sd = EmbeddedDeltaDualComplex3D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())
GLMakie.wireframe(s)

Heat = @decapode begin
    T::Form0
    D::Constant
    ∂ₜ(T) == D * ⋆(d(⋆(d(T))))
end

infer_types!(Heat, dim=3)
resolve_overloads!(Heat, dim=3)

sim = evalsim(Heat, dimension=3)

f = sim(sd, nothing, DiagonalHodge())

T_dist = MvNormal([lx/2, ly/2, lx/2], [1, 1, 1])
T = [pdf(T_dist, [p[1], p[2], p[3]]) for p in sd[:point]]

GLMakie.mesh(s; color = T, transparency=true, shading=NoShading)

D = 0.1

u₀ = ComponentArray(T=T)

tₑ = 100
prob = ODEProblem(f, u₀, (0, tₑ), (D = D,))
soln = solve(prob, Tsit5(), saveat = 1)

GLMakie.mesh(s; color = soln.u[end].T, transparency=true, shading=NoShading)

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)

  f = Figure()

  ls_T = GLMakie.LScene(f[1, 1])
  msh_T = GLMakie.mesh!(ls_T, s; color=T, colormap=:jet, transparency=true, shading=NoShading)
  Colorbar(f[1,2], msh_T)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("3d_heat.mp4", 60)