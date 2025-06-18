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
  u == ⋆(d(ψ))
  v == unit_flow(♭♯(u))

  ω == ⋆(d(u) + dᵦ(v))

  ∂ₜ(du) == -(d(⋆(∧(v, ω)))) + μ*d(⋆(d(ω)))
end

lx = ly = 8
dx = dy = 0.08
s = triangulated_grid(lx,ly,dx,dy,Point3d,true)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter());

Δ0 = Δ(0,sd);
fΔ0 = factorize(Δ0);
dual_d0 = dec_dual_derivative(0,sd);
dual_d1 = dec_dual_derivative(1,sd);
dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));

left_wall_idxs = findall(p -> p[1] < dx, s[:point]);
right_wall_idxs = findall(p -> p[1] > lx - dx, s[:point]);
bottom_wall_idxs= findall(p -> p[2] == 0, s[:point]);
top_wall_idxs = findall(p -> p[2] == ly, s[:point]);

boundary_point_idxs = union(left_wall_idxs, right_wall_idxs, bottom_wall_idxs, top_wall_idxs)

boundary_test = zeros(nv(sd))
boundary_test[boundary_point_idxs] .= 1
fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=boundary_test, colormap=:jet)
Colorbar(fig[1,2], msh)
display(fig)

top_trimmed_wall_idxs = setdiff(top_wall_idxs, union(left_wall_idxs, right_wall_idxs))

top_wall_edges = Int64[]
for e in edges(sd)
  if sd[e, :∂v0] in top_wall_idxs && sd[e, :∂v1] in top_wall_idxs
    push!(top_wall_edges, e)
  end
end

boundary_edge_test = zeros(ne(sd))
lengths = -sd[top_wall_edges, :length]
boundary_edge_test[top_wall_edges] .= lengths

# sharp_pp = ♯_mat(sd, AltPPSharp())
sharp_dd = ♯_mat(sd, LLSDDSharp());
flat_dp = ♭_mat(sd, DPPFlat())
function plot_dvf(sd, X; ls=1f0, title="Dual Vector Field")
  X♯ = sharp_dd * X
  # Makie will throw an error if the colorrange end points are not unique:
  f = Figure()
  ax = CairoMakie.Axis(f[1, 1], title=title)
  wireframe!(ax, sd, color=:gray95)
  extX = extrema(norm.(X♯))
  if (abs(extX[1] - extX[2]) > 1e-4)
    range = extX
    scatter!(ax, getindex.(sd[sd[:tri_center], :dual_point],1), getindex.(sd[sd[:tri_center], :dual_point],2), color = norm.(X♯), colorrange=range)
    Colorbar(f[1,2], limits=range)
  end
  arrows!(ax, getindex.(sd[sd[:tri_center], :dual_point],1), getindex.(sd[sd[:tri_center], :dual_point],2), getindex.(X♯,1), getindex.(X♯,2), lengthscale=ls)
  hidedecorations!(ax)
  f
end
plot_dvf(sd, boundary_edge_test; ls=0.1)

# no_flux_bc(x) = begin x[boundary_point_idxs] .= 0; return x; end
no_flux_bc(x) = return x;
unit_flow_bc(x) = begin x[top_wall_edges] .= lengths; return x; end
solve_poisson(x) = begin y = fΔ0 \ x; return y .- minimum(y); end

function generate(s, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => solve_poisson
    :dᵦ => x -> dᵦ * x
    :no_flux => no_flux_bc
    :unit_flow => unit_flow_bc
  _ => error("Unmatched operator $my_symbol")
  end
  return op
end;

sim = evalsim(incompressible_ns)
fₘ = sim(sd, generate, DiagonalHodge());

tₑ = 5.0

u₀ = ComponentArray(du = zeros(nv(sd)))

μ = 1e-3
prob = ODEProblem(fₘ, u₀, (0, tₑ), (μ=μ,))
soln = solve(prob, Tsit5(), dtmax = 0.01, saveat=tₑ/50.0, progress=true, progress_steps=1);

inv_hdg_0 = dec_inv_hodge_star(0,sd,DiagonalHodge())
fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=solve_poisson(inv_hdg_0 * soln.u[2].du), colormap=:jet)
Colorbar(fig[1,2], msh)
display(fig)

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  w = @lift(inv_hdg_0 * soln($time).du)

  f = Figure()

  ax_w = CairoMakie.Axis(f[1,1], title = @lift("Vorticity at Time $(round($time, digits=2)) with μ=$(μ)"))
  msh_w = mesh!(ax_w, s; color=w, colormap=:jet, colorrange=extrema(ω))
  Colorbar(f[1,2], msh_w)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
save_dynamics("Driven_Cavity.mp4", 30)
