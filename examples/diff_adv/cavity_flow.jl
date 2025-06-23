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
  du::DualForm2
  u::DualForm1
  ψ::Form0
  ω::Form0
  μ::Constant

  ψ == Δ⁻¹(⋆(du))
  u == no_flux(⋆(d(ψ)))
  v == flow_bc(♭♯(u))

  ω == ⋆(d(u) + dᵦ(v))

  ∂ₜ(du) == -(d(⋆(∧(v, ω)))) + μ*d(⋆(d(ω)))
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

Δ0 = Δ(0,sd);
fΔ0 = factorize(Δ0);
dual_d0 = dec_dual_derivative(0,sd);
dual_d1 = dec_dual_derivative(1,sd);
dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));

boundary_edges = findall(x -> x != 0, dual_d0 * ones(ntriangles(sd)))

function bound_edges(s, ∂₀)
  te = vcat(incident(s, ∂₀, :∂v1)...)
  se = vcat(incident(s, ∂₀, :∂v0)...)
  intersect(te, se)
end

function adj_edges(s, ∂₀)
  te = vcat(incident(s, ∂₀, :∂v1)...)
  se = vcat(incident(s, ∂₀, :∂v0)...)
  unique(vcat(te, se))
end

boundary_points = unique(vcat(s[boundary_edges, :∂v0], s[boundary_edges, :∂v1]))

boundary = zeros(nv(sd))
boundary[boundary_points] .= 2

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
CairoMakie.mesh!(ax, s; color = boundary)
fig

boundary_buffer_edges = adj_edges(sd, boundary_points)

# Dual boundary condition, boundary duals share same index as boundary primal
no_flux_bc(x) = begin x[boundary_edges] .= 0; return x; end

unit_edges = Int64[]
for e in edges(sd)
  if sd[sd[e, :∂v0], :point][2] == ly
    push!(unit_edges, e)
  end
end

lengths = sd[unit_edges, :length]

function flow_bc(x)
  x[boundary_edges] .= 0
  x[unit_edges] .= lengths
  return x
end

solve_poisson(x) = begin y = fΔ0 \ x; return y .- minimum(y); end

function generate(s, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => solve_poisson
    :dᵦ => x -> dᵦ * x
    :no_flux => no_flux_bc
    :flow_bc => flow_bc
  _ => error("Unmatched operator $my_symbol")
  end
  return op
end;

sim = evalsim(incompressible_ns)
fₘ = sim(sd, generate, GeometricHodge());

tₑ = 10.0

u₀ = ComponentArray(du = zeros(nv(sd)))

Re = 100
μ = 1/Re
prob = ODEProblem(fₘ, u₀, (0, tₑ), (μ=μ,))
soln = solve(prob, Tsit5(), dtmax = 0.01, progress=true, progress_steps=1);
soln.retcode

inv_hdg_0 = dec_inv_hodge_star(0,sd,GeometricHodge())
d0 = dec_differential(0, sd)
hdg_1 = dec_hodge_star(1, sd, GeometricHodge())
sharp_dd = ♯_mat(sd, LLSDDSharp())
flat_pd = ♭_mat(sd, DPPFlat())
sharp_pp = ♯_mat(sd, PPSharp())

function u_velocity_vector(du)
  ψ = solve_poisson(inv_hdg_0 * du)
  u = no_flux_bc(hdg_1 * (d0 * ψ))
  return sharp_dd * u
end

function v_velocity_vector(du)
  return sharp_pp * flow_bc(only.(flat_pd * u_velocity_vector(du)))
end


interp = d0_p0_interpolation(sd, hodge = GeometricHodge())
tricenter_points = sd[sd[:tri_center], :dual_point]
primal_points = sd[:point]

dual_iter = 1:20:ntriangles(sd)
primal_iter = 1:10:nv(sd)

sol = soln(3).du

u_velocity = u_velocity_vector(sol)

x = map(p -> p[1], tricenter_points[dual_iter])
y = map(p -> p[2], tricenter_points[dual_iter])
u = map(p -> p[1], u_velocity[dual_iter])
v = map(p -> p[2], u_velocity[dual_iter])

fig = Figure();
ax1 = CairoMakie.Axis(fig[1,1])
CairoMakie.arrows!(ax1, x, y, u, v; lengthscale = 0.2)

v_velocity = u_velocity_vector(sol)

x = map(p -> p[1], primal_points[primal_iter])
y = map(p -> p[2], primal_points[primal_iter])
u = map(p -> p[1], v_velocity[primal_iter])
v = map(p -> p[2], v_velocity[primal_iter])

ax2 = CairoMakie.Axis(fig[1,2])
CairoMakie.arrows!(ax2, x, y, u, v; lengthscale = 0.2)
fig


function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  psi = @lift(norm.(velocity_vector(soln($time).du)))
  velocity = @lift(map(v -> Point2d(v[1], v[2]), velocity_vector(soln($time).du)))
  u = @lift(map(p -> p[1], $velocity[iter]))
  v = @lift(map(p -> p[2], $velocity[iter]))

  f = Figure()

  ax_psi = CairoMakie.Axis(f[1,1], title = @lift("Velocity magnitude at Time $(round($time, digits=2)) with μ=$(μ)"))
  msh_psi = mesh!(ax_psi, s; color=psi, colormap=:jet)
  # CairoMakie.arrows!(ax_psi, x, y, u, v; lengthscale = 0.02)
  Colorbar(f[1,2], msh_psi)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics("Driven_Cavity_VelMag.mp4", 30)

function save_dynamics_2(save_file_name, video_length = 30)
  time = Observable(0.0)

  curl = @lift(inv_hdg_0 * soln($time).du)

  f = Figure()

  ax_curl = CairoMakie.Axis(f[1,1], title = @lift("Vorticity at Time $(round($time, digits=2)) with μ=$(μ)"))
  msh_curl = mesh!(ax_curl, s; color=curl, colormap=:jet)
  Colorbar(f[1,2], msh_curl)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics_2("Driven_Cavity_Vort.mp4", 30)
