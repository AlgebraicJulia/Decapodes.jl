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
    i += 1
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

# GMRES, Direct_LU, MG, GMRES_MG, AMG
solvername = "MG"

subs = 6
series = PrimalGeometricMapSeries(s, BinarySubdivision(), subs, Barycenter());
s = repeated_subdivision(s, binary_subdivision, subs);

md = MGData(series, sd -> Δ(0,sd), 3, BinarySubdivision());
sd = finest_mesh(series);

Δ0 = first(md.operators)

function _multigrid_μ_cycle(u, b, md::MultigridData, alg=cg, μ=1)
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
T_dist = MvNormal([lx/2.0, ly/2.0], [1/sqrt(2), 1/sqrt(2)])
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
soln = solve(prob, Tsit5(); saveat = 0.005, progress=true, progress_steps=1);

# For plotting Temperature, Pressure, Advection, and Diffusion
wdg10 = dec_wedge_product(Tuple{1, 0}, sd)
codif_1 = δ(1, sd)

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)
  f = Figure()

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
