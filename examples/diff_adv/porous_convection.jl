using Pkg
Pkg.add(["AlgebraicMultigrid",
"CairoMakie",
"Catlab",
"ComponentArrays",
"CSV",
"DataFrames",
"DiagrammaticEquations",
"Distributions",
"GeometryBasics",
"IterativeSolvers",
"Krylov",
"KrylovPreconditioners",
"LinearAlgebra",
"MAT",
"MLStyle",
"OrdinaryDiffEq",
"SparseArrays",
"StaticArrays",
"StatsBase",
"LoggingExtras",
"Logging",
"TerminalLoggers"])

using AlgebraicMultigrid
using CairoMakie
using Catlab
using CombinatorialSpaces
using ComponentArrays
using CSV
using DataFrames
using Decapodes
using DiagrammaticEquations
using Distributions
using GeometryBasics: Point2, Point3
using IterativeSolvers
using Krylov
using KrylovPreconditioners
using LinearAlgebra
using MAT
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using StaticArrays
using StatsBase
using LoggingExtras
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

import CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form
import CombinatorialSpaces.Multigrid: _multigrid_μ_cycle, car, cdr
import LinearAlgebra: ldiv!

# This model is based on the equations, constants and meshing given at
# https://pde-on-gpu.vaw.ethz.ch/lecture4

Porous_Convection = @decapode begin
  (λ_ρ₀Cp, αρ₀, k_ηf, ϕ)::Constant
  (P, T, Adv, bound_Ṫ)::Form0
  (g, qD)::Form1

  # bound_T == adiabatic(T)
  # Darcy flux
  ρ == g ∧ (αρ₀ * T)
  P == Δ⁻¹(δ(ρ))
  qD == -k_ηf * (d(P) - ρ)

  Adv == ⋆(interpolate(∧ᵈᵖ₁₁(⋆(d(T)), qD)))
  Ṫ == -1/ϕ * Adv + λ_ρ₀Cp * Δ(T)

  bound_Ṫ == tb_bc(Ṫ)

  ∂ₜ(T) == bound_Ṫ
end
infer_types!(Porous_Convection)
resolve_overloads!(Porous_Convection)
# to_graphviz(Porous_Convection)

function repeated_subdivision(s, subdivider, n)
  msh = s
  for i in 1:n
    msh = subdivider(msh)
  end
  msh
end

lx, ly = 40.0, 20.0

begin
s = EmbeddedDeltaSet2D{Bool, Point{3, Float64}}();
xs = range(0, lx; length = 3)
ys = range(0, ly; length = 3)

add_vertices!(s, 9)
i = 1
for y in ys
  for x in xs
    s[i, :point] = Point3([x, y, 0.0])
    global i += 1
  end
end
glue_sorted_triangle!(s, 1, 2, 4)
glue_sorted_triangle!(s, 2, 4, 5)
glue_sorted_triangle!(s, 2, 5, 6)
glue_sorted_triangle!(s, 2, 3, 6)
glue_sorted_triangle!(s, 4, 5, 7)
glue_sorted_triangle!(s, 5, 7, 8)
glue_sorted_triangle!(s, 5, 8, 9)
glue_sorted_triangle!(s, 5, 6, 9)

end
coarsest_domain = copy(s)

subs = 6
series = PrimalGeometricMapSeries(s, BinarySubdivision(), subs, Barycenter());
s = repeated_subdivision(s, binary_subdivision, subs);

md = MGData(series, sd -> Δ(0,sd), 3, BinarySubdivision());
sd = finest_mesh(series);

function plot_pc_domain()
  fig = Figure(size = (1600, 500),
               fontsize = 30,
               fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax1 = CairoMakie.Axis(fig[1,1],
                       title="Coarse Porous Convection Domain",
                       aspect = AxisAspect(2));
  wireframe!(ax1, coarsest_domain);
  
  ax2 = CairoMakie.Axis(fig[1,2],
                       title="Porous Convection Domain Subdivided Once",
                       aspect = AxisAspect(2),
                       yticklabelsvisible = false);
  wireframe!(ax2, binary_subdivision(coarsest_domain));
  fig
end
save("pc_coarse_domains.svg", plot_pc_domain())

# GMRES, Direct_LU, MG, GMRES_MG, AMG
solvername = "MG"
solvername = "Direct_LU"
solvername = "GMRES"

Δ0 = first(md.operators)
Δn = md.operators[begin+1]

function plot_sparsities()
  fig = Figure(size = (1700, 1200),
               fontsize = 31,
               fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax1 = CairoMakie.Axis(fig[1,1],
                       title="Sparsity Pattern of Finest Discrete Laplacian",
                       xlabel="nz = $(length(Δ0.nzval))",
                       aspect=AxisAspect(1));
  spy_plt = spy!(ax1, Δ0, markersize=10, marker=FastPixel());
  hidedecorations!(ax1, label=false)
  
  ax2 = CairoMakie.Axis(fig[1,2],
                       title="Sparsity Pattern of Second-Finest Discrete Laplacian",
                       xlabel="nz = $(length(Δn.nzval))",
                       aspect=AxisAspect(1));
  spy_plt = spy!(ax2, Δn, markersize=10, marker=FastPixel());
  hidedecorations!(ax2, label=false)
  fig
end
save("lap_sparsities.svg", plot_sparsities())

function _multigrid_μ_cycle(u, b, md::MultigridData, alg=cg, μ=1)
  A,r,p,s = CombinatorialSpaces.Multigrid.car(md)
  # Conform to CombinatorialSpaces issue #167.
  u = s == 0 ? u : alg(u,A,b,maxiter=s)
  length(md) == 1 && return u
  r_f = b - A*u
  r_c = r * r_f
  z = _multigrid_μ_cycle(zeros(size(r_c)), r_c, cdr(md), alg, μ)
  if μ > 1
    z = _multigrid_μ_cycle(z, r_c, cdr(md), alg, μ-1)
  end
  u += p * z
  # Conform to CombinatorialSpaces issue #167.
  u = s == 0 ? u : alg(u,A,b,maxiter=s)
end

function multi_solve(x; atol = √eps(Float64), rtol = √eps(Float64))
  max_iter = 50

  residuals = Float64[]
  
  y = zeros(nv(sd))
  ares(y) = norm(Δ0 * y - x)
  rNorm = ares(y)
  ϵ = atol + rtol * rNorm

  push!(residuals, rNorm)
  count = 1
  while !(ares(y) ≤ ϵ || count ≥ max_iter)
    y = _multigrid_μ_cycle(y,x,md,gauss_seidel!,2)
    push!(residuals, ares(y))
    count += 1
  end
  y, (niter = count, residuals = residuals)
end

function ldiv!(Y, A::MultigridData, B)
  Y, _ = multi_solve(B)
end  

if solvername == "Direct_LU"
  fΔ0 = LinearAlgebra.factorize(Δ0);
  Δ⁻¹(x) = begin y = fΔ0 \ x; y .-= minimum(y); end
elseif solvername == "GMRES"
  pΔ0 = ilu(Δ0);
  Δ⁻¹(x) = begin y, _ = Krylov.gmres(Δ0, x; M = pΔ0, ldiv = true, rtol = 1e-10); y .-= minimum(y) end
elseif solvername == "MG"
  Δ⁻¹(x) = begin y, _ = multi_solve(x; rtol = 1e-12, atol = 1e-10); y .-= minimum(y) end
elseif solvername == "GMRES_MG"
  Δ⁻¹(x) = begin y, _ = Krylov.gmres(Δ0, x; M = md, ldiv = true); y .-= minimum(y) end
end

mat = p2_d2_interpolation(sd)

# For no change conditions in top and bottom cooling/heating elements
bottom_wall_idxs= findall(p -> p[2] == 0, s[:point]);
top_wall_idxs = findall(p -> p[2] == ly, s[:point]);

apply_tb_bc(x) = begin x[bottom_wall_idxs] .= 0; x[top_wall_idxs] .= 0; return x; end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => Δ⁻¹
    :tb_bc => apply_tb_bc
    :interpolate => x -> mat * x
    _ => error("No operator $my_symbol found.")
  end
  return op
end

sim = eval(gensim(Porous_Convection))
f = sim(sd, generate, GeometricHodge())

ΔT = 200.0

# Initial heat distribution
# Gaussian heat disturbance along with top and bottom elements
# XXX: This constructor (with a vector as the second argument) is deprecated.
T_dist = MvNormal([lx/2.0, ly/2.0], [1/sqrt(2), 1/sqrt(2)])
# XXX: Scaling by 2 * ΔT = 400 results in a maximum temperature of ~127. It would be simpler to present were this 100.
T = [2 * ΔT * pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]
T[top_wall_idxs] .= -ΔT/2
T[bottom_wall_idxs] .= ΔT/2

# Measure the force of gravity in the downwards direction
accl_g = 9.81
grav = SVector{3}([0.0, -accl_g, 0.0])
g = eval_constant_primal_form(sd, grav)
u₀ = ComponentArray(T=T, g=g)

# Physical constants
Ra = 1000
k_ηf = 1.0
αρ₀ = (1.0/accl_g)
ϕ = 0.1
λ_ρ₀Cp = 1/Ra*(accl_g*αρ₀*k_ηf*ΔT*ly/ϕ)
constants = (k_ηf = k_ηf, αρ₀ = αρ₀, ϕ = ϕ, λ_ρ₀Cp = λ_ρ₀Cp)

tₑ = 0.5
prob = ODEProblem(f, u₀, (0, tₑ), constants)
t0 = time()
soln_decgmg = solve(prob, Tsit5();
  reltol=1e-9,
  saveat = 0.005, progress=true, progress_steps=1);
tend = time()
time_decgmg =  tend - t0

matwrite("all_solns_aug13_2025.mat",
         Dict(
              "directlu" => temperature_matrix(soln_directlu),
              "decgmg" => temperature_matrix(soln_decgmg),
              "gmres" => temperature_matrix(soln_gmres),
              "nf_directlu" => (soln_directlu.stats.nf),
              "nf_decgmg" =>   (soln_decgmg.stats.nf),
              "nf_gmres" =>    (soln_gmres.stats.nf),
              "time" => soln_directlu.t,
              "time_decgmg" => time_decgmg,
              "time_directlu" => time_directlu,
              "time_gmres" => time_gmres,
              "reltol" => 1e-9))

all_solns_mat = matread("all_solns_aug13_2025.mat")
directlu = all_solns_mat["directlu"]
decgmg = all_solns_mat["decgmg"]
gmres = all_solns_mat["gmres"]
sol_time = all_solns_mat["time"]

# For plotting Temperature, Pressure, Advection, and Diffusion
wdg10 = dec_wedge_product(Tuple{1, 0}, sd)
codif_1 = δ(1, sd)

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)
  f = Figure(size = (1600, 1200),
             fontsize = 30,
             fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));

  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Temperature at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s; color=T, colormap=:jet, colorrange=(-ΔT/2, ΔT/2))
  Colorbar(f[1,2], msh_T)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

filename = "Porous_Convection_$(solvername)_subs=$(subs)_Ra=$(Ra)"
save_dynamics("$(filename).mp4", length(soln.t))

function save_error(save_file_name, soln1, soln2, solvers::String, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln1($time).T - soln2($time).T)
  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Difference in Temp ($(solvers)) at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s; color=T, colormap=:jet, colorrange = (-0.01, 0.01))
  Colorbar(f[1,2], msh_T)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

comp1 = "LU MG"
comp2 = "LU GMRES"
comp3 = "GMRES MG"

save_error("Porous_Convection_$(comp1)_subs=$(subs)_Ra=$(Ra).mp4", soln_lu, soln_mg, comp1, length(soln.t))
save_error("Porous_Convection_$(comp2)_subs=$(subs)_Ra=$(Ra).mp4", soln_lu, soln_gmres, comp2, length(soln.t))
save_error("Porous_Convection_$(comp3)_subs=$(subs)_Ra=$(Ra).mp4", soln_gmres, soln_mg, comp3, length(soln.t))

iterations = Int64[]
for t in soln.t
  rho = codif_1 * wdg10(g, constants.αρ₀ * soln(t).T)
  _, stats = Krylov.gmres(Δ0, rho; M = pΔ0, ldiv = true, rtol = 1e-12, atol = 1e-10)
  # _, stats = multi_solve(rho; rtol = 1e-12, atol = 1e-10)
  push!(iterations, stats.niter)
end

save("$(solvername)_iterations_histo.png", hist(iterations))
save("$(solvername)_iterations_lines.png", plot(iterations))

# df = DataFrame(soln);
# CSV.write("$(filename).csv", df)

# Can use to read out rows
# out_df = CSV.read("$(filename).csv", DataFrame)
# outT = Array(out_df[2, 2:end])
#

function dynamics_composite()
  f = Figure(size = (1200, 1800),
             fontsize = 40,
             fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  ax_T0 = CairoMakie.Axis(f[1,1],
                         title = "Temperature at Time 0",
                         aspect = AxisAspect(2),
                         xticklabelsvisible = false);
  ax_T5 = CairoMakie.Axis(f[2,1],
                         title = "Temperature at Time 0.5",
                         aspect = AxisAspect(2),
                         xticklabelsvisible = false);
  ax_decgmg_diff = CairoMakie.Axis(f[3,1],
                         title = "Direct LU - DEC GMG",
                         aspect = AxisAspect(2),
                         xticklabelsvisible = false);
  ax_gmres_diff = CairoMakie.Axis(f[4,1],
                         title = "Direct LU - GMRES",
                         aspect = AxisAspect(2));

  msh_T0 = mesh!(ax_T0, s;
                 color=directlu[:,1],
                 colormap=:jet,
                 colorrange=(-ΔT/2, ΔT/2))
  msh_T5 = mesh!(ax_T5, s;
                 color=directlu[:,end],
                 colormap=:jet,
                 colorrange=(-ΔT/2, ΔT/2))

  decgmg_diff = directlu[:,end] .- decgmg[:,end]
  gmres_diff = directlu[:,end] .- gmres[:,end]
  crange_shared =
    (min(minimum(decgmg_diff), minimum(gmres_diff)),
     max(maximum(decgmg_diff), maximum(gmres_diff)))

  msh_decgmg_diff = mesh!(ax_decgmg_diff, s;
                          color=decgmg_diff,
                          colormap=:jet,
                          colorrange=crange_shared)
  msh_gmres_diff = mesh!(ax_gmres_diff, s;
                         color=gmres_diff,
                         colormap=:jet,
                         colorrange=crange_shared)

  Colorbar(f[1:2,2], msh_T0, height=Relative(1.0))
  Colorbar(f[3:4,2], msh_decgmg_diff, height=Relative(1.0))
  f
end
dynamics_composite()

filename = "Porous_Convection_$(solvername)_subs=$(subs)_Ra=$(Ra)"
save_dynamics("$(filename).mp4", length(soln.t))

function save_error(save_file_name, soln1, soln2, solvers::String, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln1($time).T - soln2($time).T)
  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Difference in Temp ($(solvers)) at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s; color=T, colormap=:jet, colorrange = (-0.01, 0.01))
  Colorbar(f[1,2], msh_T)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

comp1 = "LU MG"
comp2 = "LU GMRES"
comp3 = "GMRES MG"

save_error("Porous_Convection_$(comp1)_subs=$(subs)_Ra=$(Ra).mp4", soln_lu, soln_mg, comp1, length(soln.t))
save_error("Porous_Convection_$(comp2)_subs=$(subs)_Ra=$(Ra).mp4", soln_lu, soln_gmres, comp2, length(soln.t))
save_error("Porous_Convection_$(comp3)_subs=$(subs)_Ra=$(Ra).mp4", soln_gmres, soln_mg, comp3, length(soln.t))

iterations = Int64[]
for t in soln.t
  rho = codif_1 * wdg10(g, constants.αρ₀ * soln(t).T)
  _, stats = Krylov.gmres(Δ0, rho; M = pΔ0, ldiv = true, rtol = 1e-12, atol = 1e-10)
  # _, stats = multi_solve(rho; rtol = 1e-12, atol = 1e-10)
  push!(iterations, stats.niter)
end

save("$(solvername)_iterations_histo.png", hist(iterations))
save("$(solvername)_iterations_lines.png", plot(iterations))

# df = DataFrame(soln);
# CSV.write("$(filename).csv", df)

# Can use to read out rows
# out_df = CSV.read("$(filename).csv", DataFrame)
# outT = Array(out_df[2, 2:end])

# SVGs of such plots will exhibit fuzzy text.

function dynamics_FCS(soln, titleextra)
  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1],
                         title = "Temperature at Time 0.5 ; $titleextra",
                         aspect = AxisAspect(2))
  msh_T = mesh!(ax_T, s;
                color=soln(0.5).T,
                colormap=:jet,
                colorrange=(-ΔT/2, ΔT/2))

  # https://discourse.julialang.org/t/makie-jl-limiting-height-of-colorbar-to-that-of-a-nearby-equal-aspect-axis/55876/4
  # Pixel-accurate value is .6351:
  Colorbar(f[1,2], msh_T, height=Relative(0.6351))

  f
end
save("temperature_ics.png", dynamics_ICS(soln_decgmg))

function dynamics_FCS(soln, titleextra)
  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1],
                         title = "Temperature at Time 0.5 ; $titleextra",
                         aspect = AxisAspect(2))
  msh_T = mesh!(ax_T, s;
                color=soln(0.5).T,
                colormap=:jet,
                colorrange=(-ΔT/2, ΔT/2))

  # https://discourse.julialang.org/t/makie-jl-limiting-height-of-colorbar-to-that-of-a-nearby-equal-aspect-axis/55876/4
  # Pixel-accurate value is .6351:
  Colorbar(f[1,2], msh_T, height=Relative(0.6351))

  f
end
dynamics_FCS(soln_decgmg, "DEC GMG")

# SVGs of such plots will exhibit fuzzy text.
save("temperature_fcs_decgmg.png", dynamics_FCS(soln_decgmg, "DEC GMG"))

function diffs_FCS(soln, titleextra)
  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1],
                         title = "Temperature Difference at Time 0.5 ; Direct LU - $titleextra",
                         aspect = AxisAspect(2))

  msh_T = mesh!(ax_T, s;
                color=soln_directlu[:,end] - soln[:,end],
                colormap=:jet)

  # https://discourse.julialang.org/t/makie-jl-limiting-height-of-colorbar-to-that-of-a-nearby-equal-aspect-axis/55876/4
  # Pixel-accurate value is .605:
  Colorbar(f[1,2], msh_T, height=Relative(0.605))

  f
end
diffs_FCS(soln_gmres, "GMRES")

# SVGs of such plots will exhibit fuzzy text.
save("temperature_diff_decgmg.png", diffs_FCS(soln_decgmg, "DEC GMG"))
save("temperature_diff_gmres.png", diffs_FCS(soln_gmres, "GMRES"))

sqrt(mean((soln_directlu(0.5).T - soln_gmres(0.5).T) .^ 2))

rmse(ana_sol, num_sol) = sqrt(mean((ana_sol .- num_sol) .^ 2))

function plot_rmse_over_time(num_sol, titleextra)
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="RMSE of $titleextra over Time",
                       xlabel="time [s]",
                       ylabel="RMSE [°C]")
  rmse_over_time = map(soln_directlu.t) do t
    rmse(soln_directlu(t).T, num_sol(t).T)
  end
  lines!(ax, soln_directlu.t, rmse_over_time);
  fig
end

plot_rmse_over_time(soln_decgmg, "DEC GMG")
plot_rmse_over_time(soln_gmres, "GMRES")

function plot_comparison_rmse_over_time()
  fig = Figure(size = (1600, 800),
               fontsize = 30,
               fonts = (; regular = "CMU Serif Roman", bold = "CMU Serif Bold"));
  axlogit = CairoMakie.Axis(fig[1,2],
                       title="log10 RMSE over Time",
                       xlabel="time [s]",
                       ylabel= "log10 RMSE",
                       yaxisposition = :right)
  axnologit = CairoMakie.Axis(fig[1,1],
                       title="RMSE over Time",
                       xlabel="time [s]",
                       ylabel="RMSE [°C]")

  rmse_over_time_decgmg = map(soln_directlu.t) do t
    rmse(soln_directlu(t).T, soln_decgmg(t).T)
  end
  rmse_over_time_gmres = map(soln_directlu.t) do t
    rmse(soln_directlu(t).T, soln_gmres(t).T)
  end

  lines!(axlogit, soln_directlu.t[2:end], log10.(rmse_over_time_gmres)[2:end], label="GMRES");
  lines!(axlogit, soln_directlu.t[2:end], log10.(rmse_over_time_decgmg)[2:end], label="DEC GMG");
  lines!(axnologit, soln_directlu.t, rmse_over_time_gmres, label="GMRES");
  lines!(axnologit, soln_directlu.t, rmse_over_time_decgmg, label="DEC GMG");

  #Legend(fig[1,2], ax);
  axislegend(axnologit,"Scheme", position=:lt, labelsize=24, fontsize=30);
  fig
end
plot_comparison_rmse_over_time()
save("rmse_comparison.svg", plot_comparison_rmse_over_time())

rmse_ratios = rmse_over_time_decgmg ./ rmse_over_time_gmres

function temperature_matrix(x)
  reduce(hcat, map(x) do v
    v.T
  end)
end

matwrite("rmses.mat", Dict("rmse_over_time_decgmg" => rmse_over_time_decgmg,
                           "rmse_over_time_gmres" => rmse_over_time_gmres))

using CairoMakie
using MAT

data = matread("rmses.mat")
rmse_over_time_decgmg = data["rmse_over_time_decgmg"]
rmse_over_time_gmres = data["rmse_over_time_gmres"]

f = Figure();
ax = CairoMakie.Axis(f[1,1], title="RMSE RATIOS");
ax2 = CairoMakie.Axis(f[2,1], title="(log10 RMSE) RATIOS");
lines!(ax, (rmse_over_time_decgmg ./ rmse_over_time_gmres)[2:end])
lines!(ax2, (log10.(rmse_over_time_decgmg)[2:end] ./ log10.(rmse_over_time_gmres)[2:end]));
f
save("rmse_ratios_and_ratios_of_logs.png", ans)


