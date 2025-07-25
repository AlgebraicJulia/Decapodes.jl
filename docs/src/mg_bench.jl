using ACSets
using CombinatorialSpaces
using BenchmarkTools
using LinearAlgebra
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
f = lap * rand(nv(sd));
#using MAT
#matwrite("lap_times_random.mat", Dict("f" => f))
#f = matread("lap_times_random.mat")["f"]

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
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="L r",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=f,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
plot_ics()
save("laplacian_times_random.png", plot_ics())

function plot_directlu_solve()
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="L \\ b via Direct LU solve",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
plot_directlu_solve()
save("poisson_solve_directlu.png", plot_directlu_solve())

function plot_mg_solve()
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="L \\ b via DEC GMG solve",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=mg_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
plot_mg_solve()
save("poisson_solve_mg.png", plot_mg_solve())

function plot_direct_minus_mg()
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="Direct LU - DEC GMG",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u .- mg_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
plot_direct_minus_mg()
save("poisson_directlu_minus_mg.png", plot_direct_minus_mg())

function plot_direct_minus_ilu0cg()
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="Direct LU - ILU0 CG",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u .- kry_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
plot_direct_minus_ilu0cg()
save("poisson_directlu_minus_ilu0cg.png", plot_direct_minus_ilu0cg())

function plot_direct_minus_gmgcg()
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="Direct LU - GMG CG",
                       aspect = AxisAspect(1));
  msh = mesh!(ax, sp3;
              color=ds_u .- kr_mg_u,
              colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
plot_direct_minus_gmgcg()
save("poisson_directlu_minus_gmgcg.png", plot_direct_minus_gmgcg())

