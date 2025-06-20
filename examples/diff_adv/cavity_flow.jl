using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using CairoMakie
using ComponentArrays
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays

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
  # u == ⋆(d(ψ))
  v == ♭(no_flow(unit_flow(♯(u))))

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

off_boundary_edges = findall(x -> x != 0, dual_d0 * ones(ntriangles(sd)))

no_flux_bc(x) = begin x[off_boundary_edges] .= 0; return x; end

function unit_flow_bc(x)
  nx = length(0:dx:lx)
  nt = ntriangles(sd)
  nlast_row = 2*(nx-1)
  for i in nt-(nlast_row)+1:nt
    x[i] = Point3d(1,0,0)
  end
  return x
end

function no_flow_bc(x)
  nx = length(0:dx:lx)
  ny = length(0:dy:ly)

  velocity = Point3d(0,0,0)

  trow = 2*(nx-1)
  for i in 1:trow
    x[i] = velocity
  end
  for y in 1:ny-1
    base_idx = (y-1)*trow
    x[base_idx + 1] = velocity
    x[base_idx + trow] = velocity
  end
  return x
end

X♯ = zeros(Point3d, ntriangles(sd))
unit_flow_bc(X♯)
flat_pd = ♭_mat(sd, DPPFlat())
sharp_pp = ♯_mat(sd, PPSharp())

Y♯ = sharp_pp * flat_pd * X♯

# no_flow_bc(X♯)

# function plot_vf(sd, X♯; ls=1f0, title="Dual Vector Field")
#   # Makie will throw an error if the colorrange end points are not unique:
#   f = Figure()
#   ax = CairoMakie.Axis(f[1, 1], title=title)
#   wireframe!(ax, sd, color=:gray95)
#   extX = extrema(norm.(X♯))
#   if (abs(extX[1] - extX[2]) > 1e-4)
#     range = extX
#     scatter!(ax, getindex.(sd[sd[:tri_center], :dual_point],1), getindex.(sd[sd[:tri_center], :dual_point],2), color = norm.(X♯), colorrange=range)
#     Colorbar(f[1,2], limits=range)
#   end
#   arrows!(ax, getindex.(sd[sd[:tri_center], :dual_point],1), getindex.(sd[sd[:tri_center], :dual_point],2), getindex.(X♯,1), getindex.(X♯,2), lengthscale=ls)
#   hidedecorations!(ax)
#   f
# end
# plot_vf(sd, X♯; ls=0.1)

solve_poisson(x) = begin y = fΔ0 \ x; return y .- minimum(y); end

function generate(s, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => solve_poisson
    :dᵦ => x -> dᵦ * x
    :no_flux => no_flux_bc
    :no_flow => no_flow_bc
    :unit_flow => unit_flow_bc
  _ => error("Unmatched operator $my_symbol")
  end
  return op
end;

sim = evalsim(incompressible_ns)
fₘ = sim(sd, generate, GeometricHodge());

tₑ = 3

u₀ = ComponentArray(du = zeros(nv(sd)))

Re = 100
μ = 1/Re
prob = ODEProblem(fₘ, u₀, (0, tₑ), (μ=μ,))
soln = solve(prob, Tsit5(), dtmax = 0.01, saveat = 0.01, progress=true, progress_steps=1);

inv_hdg_0 = dec_inv_hodge_star(0,sd,GeometricHodge())
d0 = dec_differential(0, sd)
hdg_1 = dec_hodge_star(1, sd, GeometricHodge())
sharp_dd = ♯_mat(sd, LLSDDSharp())
function velocity_vector(du)
  ψ = solve_poisson(inv_hdg_0 * du)
  u = no_flux_bc(hdg_1 * (d0 * ψ))
  return no_flow_bc(unit_flow_bc(sharp_dd * u))
end

velocity = map(v -> Point2d(v[1], v[2]), velocity_vector(soln.u[end].du))
tricenter_points = sd[sd[:tri_center], :dual_point]

iter = 1:10:ntriangles(sd)

x = map(p -> p[1], tricenter_points[iter])
y = map(p -> p[2], tricenter_points[iter])

u = map(p -> p[1], velocity[iter])
v = map(p -> p[2], velocity[iter])

interp = d0_p0_interpolation(sd, hodge = GeometricHodge())

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=interp * norm.(velocity_vector(soln.u[end].du)), colormap=:jet)
CairoMakie.arrows!(ax, x, y, u, v; lengthscale = 0.02)
Colorbar(fig[1,2], msh)
display(fig)

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  psi = @lift(interp * norm.(velocity_vector(soln($time).du)))

  f = Figure()

  ax_psi = CairoMakie.Axis(f[1,1], title = @lift("Velocity magnitude at Time $(round($time, digits=2)) with μ=$(μ)"))
  msh_psi = mesh!(ax_psi, s; color=psi, colormap=:jet)
  Colorbar(f[1,2], msh_psi)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics("Driven_Cavity.mp4", 30)
