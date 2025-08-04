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
ly = lz = 30

s = parallelepiped(lx = lx, ly = ly, lz = lz; tetcmd = "vpVq2a5")
sd = EmbeddedDeltaDualComplex3D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

# GLMakie.wireframe(s)

Adv = @decapode begin
    T::Form0
    U::Form1
    c::Constant
    ∂ₜ(T) == -c*(⋆(d(⋆(T ∧ U))))
end

infer_types!(Adv, dim=3)
resolve_overloads!(Adv, dim=3)

sim = evalsim(Adv, dimension=3)

f = sim(sd, nothing, DiagonalHodge())

T_dist = MvNormal([lx/2, ly/2, lz/2], [5,5,5])
T = [pdf(T_dist, [p[1], p[2], p[3]]) for p in sd[:point]]
T₀ = T
GLMakie.mesh(s; color = T, transparency=true, shading=NoShading)

U = eval_constant_primal_form(sd, Point3d(1,0,0))

u₀ = ComponentArray(T=T, U=U)

c = 1

tₑ = 20
prob = ODEProblem(f, u₀, (0, tₑ), (c=c,))
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