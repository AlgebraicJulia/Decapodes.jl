using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using CairoMakie
using ComponentArrays
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using ACSets

using LoggingExtras
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

# Reproduced from Mohamed et al. (2016), simulation 4.5
incompressible_ns = @decapode begin
  ψ::Form0
  u::DualForm1
  v::Form1
  ω::Form0
  μ::Constant

  u == ⋆(d(flow_bc(ψ)))
  v == ♭♯(u)

  ω == ⋆(d(u) + dᵦ(v))

  ∂ₜ(ψ) == no_change_bc(dsdinv(-(d(⋆(∧(v, ω)))) + μ*d(⋆(d(ω)))))
end

lx = ly = 1
dx = dy = 0.01
s = triangulated_grid(lx,ly,dx,dy,Point3d,false)
nx = length(0:dx:lx)
ny = length(0:dy:ly)

for y in 1:ny
  if iseven(y)
    for x in 1:nx
      s[(y - 1) * nx + x, :point] = s[(y - 1) * nx + x, :point] .- Point3d(dx/2,0,0)
    end
  end
end

sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Barycenter());

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
wireframe!(ax, s)
fig

dual_d0 = dec_dual_derivative(0,sd);
dual_d1 = dec_dual_derivative(1,sd);
dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));

hdg_1 = dec_hodge_star(1, sd, DiagonalHodge())
d0 = dec_differential(0, sd)
dsd = dual_d1 * hdg_1 * d0
dsdinv = factorize(dsd)

boundary_edges = findall(x -> x != 0, dual_d0 * ones(ntriangles(sd)))
boundary_points = unique(vcat(s[boundary_edges, :∂v0], s[boundary_edges, :∂v1]))

boundary_buffer_points = deepcopy(boundary_points)

for e in edges(sd)
  if sd[e, :∂v0] in boundary_points
    push!(boundary_buffer_points, sd[e, :∂v1])
  elseif sd[e, :∂v1] in boundary_points
    push!(boundary_buffer_points, sd[e, :∂v0])
  end
end

boundary_buffer_points = unique(boundary_buffer_points)

nx = length(0:dx:lx)
top_stream_points = collect(nv(sd)-2*nx+2:nv(sd)-nx-1)

# boundary = zeros(nv(sd))
# boundary[boundary_buffer_points] .= 1
# boundary[boundary_points] .= 2
# boundary[top_stream_points] .= 3

# fig = Figure();
# ax = CairoMakie.Axis(fig[1,1])
# CairoMakie.mesh!(ax, s; color = boundary)
# fig

function flow_bc(x)
  x[boundary_buffer_points] .= 0;
  x[top_stream_points] .= dy;
  return x;
end

function no_change_bc(x)
  x[boundary_buffer_points] .= 0;
  return x;
end

# boundary = flow_bc(zeros(nv(sd)))
# fig = Figure();
# ax = CairoMakie.Axis(fig[1,1])
# CairoMakie.mesh!(ax, s; color = boundary)
# fig

# solve_poisson(x) = begin y = fΔ0 \ x; return y .- minimum(y); end
solve_dsdinv(x) = begin x -> dsdinv \ x; x .= x .- x[1]; return x; end

function generate(s, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => solve_poisson
    :dᵦ => x -> dᵦ * x
    :dsdinv => solve_dsdinv
    :flow_bc => flow_bc
    :no_change_bc => no_change_bc
  _ => error("Unmatched operator $my_symbol")
  end
  return op
end;

sim = evalsim(incompressible_ns)
fₘ = sim(sd, generate, GeometricHodge());

tₑ = 10.0

u₀ = ComponentArray(ψ = zeros(nv(sd)))

Re = 100
μ = 1/Re
prob = ODEProblem(fₘ, u₀, (0, tₑ), (μ=μ,))
soln = solve(prob, Tsit5(), dtmax = 0.01, progress=true, progress_steps=1);
soln.retcode

fig = Figure();
ax_psi = CairoMakie.Axis(fig[1,1])
msh_psi = mesh!(ax_psi, s; color=soln.u[1].ψ, colormap=:jet)
Colorbar(fig[1,2], msh_psi)
fig


hdg_1 = dec_hodge_star(1, sd, GeometricHodge())
sharp_dd = ♯_mat(sd, LLSDDSharp())
# flat_pd = ♭_mat(sd, DPPFlat())
# sharp_pp = ♯_mat(sd, PPSharp())

function u_velocity_vector(psi)
  u = hdg_1 * (d0 * psi)
  return sharp_dd * u
end

function v_velocity_vector(du)
  return sharp_pp * flow_bc(only.(flat_pd * u_velocity_vector(du)))
end

interp = d0_p0_interpolation(sd, hodge = GeometricHodge())
tricenter_points = sd[sd[:tri_center], :dual_point]

dual_iter = 1:20:ntriangles(sd)

x = map(p -> p[1], tricenter_points[dual_iter])
y = map(p -> p[2], tricenter_points[dual_iter])

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  psi = @lift(soln($time).ψ)
  # velocity = @lift(map(v -> Point2d(v[1], v[2]), u_velocity_vector($psi)))
  # u = @lift(map(p -> p[1], $velocity[dual_iter]))
  # v = @lift(map(p -> p[2], $velocity[dual_iter]))

  f = Figure()

  ax_psi = CairoMakie.Axis(f[1,1], title = @lift("Streamfunction at Time $(round($time, digits=2)) with μ=$(μ)"))
  msh_psi = mesh!(ax_psi, s; color=psi, colormap=:jet)
  # CairoMakie.arrows!(ax_psi, x, y, u, v; lengthscale = 0.02)
  Colorbar(f[1,2], msh_psi)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics("Driven_Cavity_Streamfunction.mp4", 30)

# function save_dynamics_2(save_file_name, video_length = 30)
#   time = Observable(0.0)

#   curl = @lift(inv_hdg_0 * soln($time).du)

#   f = Figure()

#   ax_curl = CairoMakie.Axis(f[1,1], title = @lift("Vorticity at Time $(round($time, digits=2)) with μ=$(μ)"))
#   msh_curl = mesh!(ax_curl, s; color=curl, colormap=:jet)
#   Colorbar(f[1,2], msh_curl)

#   timestamps = range(0, soln.t[end], length=video_length)
#   record(f, save_file_name, timestamps; framerate = 15) do t
#     time[] = t
#   end
# end
# save_dynamics_2("Driven_Cavity_Vort.mp4", 30)
