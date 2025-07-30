using ACSets
using CombinatorialSpaces
using BenchmarkTools
using LinearAlgebra
using StaticArrays
using Krylov, KrylovPreconditioners
using IterativeSolvers: gauss_seidel!, ssor!
import CombinatorialSpaces.Multigrid: _multigrid_μ_cycle, multigrid_μ_cycles, multigrid_wcycles, car, cdr

function _multigrid_μ_cycle(u, b, md::MultigridData, alg=gauss_seidel!, μ=1)
  A,r,p,s = CombinatorialSpaces.Multigrid.car(md)
  u = alg(u,A,b,maxiter=s)
  length(md) == 1 && return u
  r_f = b - A*u
  r_c = r * r_f
  z = _multigrid_μ_cycle(zeros(size(r_c)), r_c, cdr(md), alg, μ)
  if μ > 1
    z = _multigrid_μ_cycle(z, r_c, cdr(md), alg, μ-1)
  end
  u += p * z
  u = alg(u,A,b,maxiter=s)
end

function multigrid_μ_cycles(u, b, md::MultigridData, cycles, alg=cg, μ=1)
  cycles == 0 && return u
  u = _multigrid_μ_cycle(u,b,md,alg,μ)
  multigrid_μ_cycles(u,b,md,cycles-1,alg,μ)
end

multigrid_wcycles(u, b, md, cycles, alg=gauss_seidel!) = multigrid_μ_cycles(u, b, md, cycles, alg, 2)

# Edit these constants for benchmarking!
LEVELS = 7
REMOVE_LAYERS = 3

println("Generating Mesh")
s = triangulated_grid(240, 240, 10, 10, Point2D);
multi_mesh = PrimalGeometricMapSeries(s, binary_subdivision_map, LEVELS - 1);
sd = finest_mesh(multi_mesh);

println("Mesh elements:")
println(" Vertices: $(nv(sd))")
println(" Edges: $(ne(sd))")
println(" Triangles: $(ntriangles(sd))")

println(" Dual vertices: $(nparts(sd, :DualV))")
println(" Dual edges: $(nparts(sd, :DualE))")
println(" Dual triangles: $(nparts(sd, :DualTri))")

lap = ∇²(0, sd)

# The initial conditions for temperature from the Porous Convection problem.
ΔT = 200.0
using Distributions
T_dist = MvNormal([120, 120], [0.5 0 ; 0 0.5])
T = [2 * ΔT * pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]
bottom_wall_idxs= findall(p -> p[2] == 0, s[:point]);
top_wall_idxs = findall(p -> p[2] == 240, s[:point]);
T[top_wall_idxs] .= -ΔT/2
T[bottom_wall_idxs] .= ΔT/2

# Gravity:
accl_g = 9.81
grav = SVector{3}([0.0, -accl_g, 0.0])
g = eval_constant_primal_form(sd, grav)
# αρ₀:
αρ₀ = (1.0/accl_g)

wdg10 = dec_wedge_product(Tuple{1,0}, sd)
s1 = dec_hodge_star(1, sd, GeometricHodge())
dd_1 = dec_dual_derivative(1, sd)
is0 = dec_inv_hodge_star(0,sd)
rho = wdg10(g, αρ₀ .* T)
# This is the initial δ(ρ) as it appears in the Porous Convection Decapode.
div_rho = is0 * dd_1 * s1 * rho

f = div_rho

#using MAT
#matwrite("div_rho.mat", Dict("div_rho" => div_rho))

#fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
#ax = CairoMakie.Axis(fig[1,1]);
#mesh!(ax, sp3, color=(log10 ∘ abs).(div_rho));
#fig

#f = lap * rand(nv(sd));
#using MAT
#matwrite("lap_times_random.mat", Dict("f" => f))
#f = matread("lap_times_random.mat")["f"]

#f = matread("div_rho.mat")["div_rho"]

# Steady state Δu + f = 0 => Δu = -f => u = -∇⁻¹(f)

rel_err(A, x, b; p = 2) = norm(b - A * x, p) / norm(b, p)

# Direct Solver
println("Direct Solver (LU)")
lap_lu = lu(lap)

ds_u = lap_lu \ f;
@btime ds_u = lap_lu \ f;
ds_u = lap_lu \ f;
@show rel_err(lap, ds_u, f)

# Krylov with ilu(0)
println("ILU0 Preconditioned CG")
lap_ilu0 = ilu(lap)

@btime kry_u, stats = cg(lap, f; M = lap_ilu0, ldiv = true);
kry_u, stats = cg(lap, f; M = lap_ilu0, ldiv = true);
@show rel_err(lap, kry_u, f)

cdr2(md::MultigridData) =
  length(md) > 1 ?
    MultigridData(md.operators[1:end-1],md.restrictions[1:end-1],md.prolongations[1:end-1],md.steps[1:end-1]) :
    error("Not enough grids remaining in $md to take the cdr2.")

import LinearAlgebra: ldiv!

function ldiv!(y, A::MultigridData, b)
  y .= multigrid_wcycles(zeros(nv(sd)), b, A, 2);
end

md = MGData(multi_mesh, sd -> ∇²(0, sd), 3, BinarySubdivision());
for i in 1:REMOVE_LAYERS+1
  # CombinatorialSpaces Multigrid
  println("Multigrid Solver, with $(LEVELS - i + 1) layers")
  @btime mg_u = multigrid_wcycles(zeros(nv(sd)), f, md, 5);
  mg_u = multigrid_wcycles(zeros(nv(sd)), f, md, 5);
  @show rel_err(lap, mg_u, f)

  println("Geometric Multigrid Preconditioned CG, with $(LEVELS - i + 1) layers")

  @btime kr_mg_u, kr_mg_stats = cg(lap, f, M = md, ldiv = true);
  kr_mg_u, kr_mg_stats = cg(lap, f, M = md, ldiv = true);
  @show rel_err(lap, kr_mg_u, f)

  global md = cdr2(md)
end

sp3 = triangulated_grid(240, 240, 10, 10, Point3D);
sp3 = binary_subdivision(binary_subdivision(binary_subdivision(binary_subdivision(binary_subdivision(binary_subdivision(sp3))))))

using CairoMakie
function plot_ics()
  fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax = CairoMakie.Axis(fig[1,1],
                       title="δ(ρ)",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=f,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
#plot_ics()
save("divrho.png", plot_ics())

function plot_directlu_solve()
  fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax = CairoMakie.Axis(fig[1,1],
                       title="L \\ b via Direct LU solve",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
#plot_directlu_solve()
save("poisson_solve_directlu_divrho.png", plot_directlu_solve())

function plot_mg_solve()
  fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax = CairoMakie.Axis(fig[1,1],
                       title="L \\ b via DEC GMG solve",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=mg_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
#plot_mg_solve()
save("poisson_solve_mg_divrho.png", plot_mg_solve())

function plot_direct_minus_mg()
  fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax = CairoMakie.Axis(fig[1,1],
                       title="Direct LU - DEC GMG",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=(ds_u .- mean(ds_u)) .- (mg_u .- mean(mg_u)),
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
#plot_direct_minus_mg()
save("poisson_directlu_minus_mg_divrho.png", plot_direct_minus_mg())

function plot_direct_minus_ilu0cg()
  fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax = CairoMakie.Axis(fig[1,1],
                       title="Direct LU - ILU0 CG",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u .- kry_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
#plot_direct_minus_ilu0cg()
save("poisson_directlu_minus_ilu0cg_divrho.png", plot_direct_minus_ilu0cg())

function plot_direct_minus_gmgcg()
  fig = Figure(fontsize = 20, fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax = CairoMakie.Axis(fig[1,1],
                       title="Direct LU - GMG CG",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u .- kr_mg_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
#plot_direct_minus_gmgcg()
save("poisson_directlu_minus_gmgcg_divrho.png", plot_direct_minus_gmgcg())

function plot_main_composite(centered=true)
  fig = Figure(size = (1600, 1200),
               fontsize = 30,
               fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax_topleft = CairoMakie.Axis(fig[1,2],
                       title = "L \\ b via Direct LU solve",
                       aspect = AxisAspect(1),
                       xticklabelsvisible = false);
  ax_toprght = CairoMakie.Axis(fig[1,4],
                       title = "Direct LU - DEC GMG",
                       aspect = AxisAspect(1),
                       xticklabelsvisible = false,
                       yticklabelsvisible = false);
  ax_botleft = CairoMakie.Axis(fig[2,2],
                       title = "Direct LU - GMG CG",
                       aspect = AxisAspect(1));
  ax_botrght = CairoMakie.Axis(fig[2,4],
                       title = "Direct LU - ILU0 CG",
                       aspect = AxisAspect(1),
                       yticklabelsvisible = false);

  direct_msh = mesh!(ax_topleft, sp3;
              color=ds_u,
              colormap=:jet)

  centered_diff(x) = (ds_u .- mean(ds_u)) .- (x .- mean(x))
  difference = centered ? centered_diff : (x -> ds_u .- x)
  decgmg_diffs = difference(mg_u)
  gmgcg_diffs  = difference(kr_mg_u)
  kry_u_diffs  = difference(kry_u)

  crange_shared =
    (min(minimum(centered_gmgcg_diffs), minimum(centered_kry_u_diffs)),
     max(maximum(centered_gmgcg_diffs), maximum(centered_kry_u_diffs)))

  min_decgmg_msh = mesh!(ax_toprght, sp3;
              color=decgmg_diffs,
              colormap=:jet)

  min_gmgcg_msh = mesh!(ax_botleft, sp3;
              color=gmgcg_diffs,
              colormap=:jet,
              colorrange = crange_shared)

  min_ilu0cg_msh = mesh!(ax_botrght, sp3;
              color=kry_u_diffs,
              colormap=:jet,
              colorrange = crange_shared)

  Colorbar(fig[1,1], direct_msh)
  Colorbar(fig[1,5], min_decgmg_msh)
  Colorbar(fig[2,3], min_gmgcg_msh)
  fig
end
plot_main_composite()
save("pc_main_composite.png", plot_main_composite())

using MAT
matwrite("divrho_solutions.mat", Dict([
  "lap" => lap,
  "f" => div_rho,
  "ds_u" => ds_u,
  "mg_u" => mg_u,
  "kr_mg_u" => kr_mg_u,
  "kry_u" => kry_u]))

