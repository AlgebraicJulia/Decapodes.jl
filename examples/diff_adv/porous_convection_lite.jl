#=
using Pkg
Pkg.add(["AlgebraicMultigrid",
"CairoMakie",
"Catlab",
"CombinatorialSpaces",
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
=#

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

lx, ly = 40.0, 20.0

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

subs = 6
series = PrimalGeometricMapSeries(s, BinarySubdivision(), subs, Barycenter());

md = MGData(series, sd -> Δ(0,sd), 3, BinarySubdivision());
sd = finest_mesh(series);

Δ0 = first(md.operators)

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
bottom_wall_idxs= findall(p -> p[2] == 0, sd[:point]);
top_wall_idxs = findall(p -> p[2] == ly, sd[:point]);

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
soln = solve(prob, Tsit5();
             saveat = 0.005, progress=true, progress_steps=1);

@show (soln.retcode, soln.stats.nf, soln.stats.naccept, soln.stats.nreject)

function temperature_matrix(x)
  reduce(hcat, map(x) do v
    v.T
  end)
end

matwrite("soln.mat",
  Dict(["soln"    => temperature_matrix(soln),
        "time"    => soln.t,
        "nf"      => soln.stats.nf,
        "naccept" => soln.stats.naccept,
        "nreject" => soln.stats.nreject])
soln = matread("soln.mat")["soln"]

all_solns_mat = matread("all_solns.mat")
directlu = all_solns_mat["directlu"]
sol_time = all_solns_mat["time"]

function plot_rmse_over_time()
  fig = Figure();
  ax = CairoMakie.Axis(fig[1,1],
                       title="RMSE over Time",
                       xlabel="time [s]",
                       ylabel="RMSE [°C]")
  rmse_over_time = map(eachindex(sol_time)) do ti
    rmse(directlu[:,ti], soln[:,ti])
  end
  @show rmse_over_time[end]
  lines!(ax, sol_time, rmse_over_time);
  fig
end

save("RMSE_over_time.svg", plot_rmse_over_time(soln)

