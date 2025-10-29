using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using ComponentArrays
using OrdinaryDiffEq
using GLMakie
using Distributions
using SparseArrays

lx = 100
ly = lz = 30

s = parallelepiped(lx = lx, ly = ly, lz = lz; tetcmd = "vpVq2a5")
sd = EmbeddedDeltaDualComplex3D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Barycenter())

# GLMakie.wireframe(s)

# Adv = @decapode begin
#     T::Form0
#     U::Form1
#     c::Constant
#     ∂ₜ(T) == -c*(δ(T ∧ U))
# end

Adv = @decapode begin
    T::DualForm0
    U::Form1
    c::Constant
    ∂ₜ(T) == -c*⋆(∧(⋆(d(T)), U))
end

infer_types!(Adv, dim=3)
resolve_overloads!(Adv, dim=3)

gensim(Adv, dimension=3)

sim = evalsim(Adv, dimension=3)

f = sim(sd, nothing, DiagonalHodge())

T_dist = MvNormal([lx/2, ly/2, lz/2], [5,5,5])
# T = [pdf(T_dist, [p[1], p[2], p[3]]) for p in sd[:point]]
T = [pdf(T_dist, [p[1], p[2], p[3]]) for p in sd[tetrahedron_center(sd), :dual_point]]

T₀ = T

function p3_d3_interpolation(sd::HasDeltaSet3D)
  mat = spzeros(nv(sd), ntetrahedra(sd))
  for tet_id in tetrahedra(sd)
    tet_vol = sd[tet_id, :vol]
    for dual_tet_id in (1:24) .+ 24 * (tet_id - 1)
      dual_tet_vol = sd[dual_tet_id, :dual_vol]

      weight = (dual_tet_vol / tet_vol)

      v = sd[sd[sd[dual_tet_id, :D_∂t1], :D_∂e2], :D_∂v1]

      mat[v, tet_id] += weight
    end
  end

  mat
end

interp = d0_p0_interpolation(sd)

GLMakie.mesh(s; color = interp * T, transparency=true, shading=NoShading)

U = eval_constant_primal_form(sd, Point3d(1,0,0))

u₀ = ComponentArray(T=T, U=U)

c = 1

tₑ = 20
prob = ODEProblem(f, u₀, (0, tₑ), (c=c,))
soln = solve(prob, Tsit5(), saveat = 0.1)

GLMakie.mesh(s; color = interp * soln.u[end].T, transparency=true, shading=NoShading, colorrange = extrema(interp * T₀))

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