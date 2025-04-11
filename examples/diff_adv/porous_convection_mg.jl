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

function repeated_subdivision(s, subdivider, n)
  msh = s
  for i in 1:n
    msh = subdivider(msh)
  end
  msh
end

lx, ly = 22.0, 20.0
dx = dy = 20

s = triangulated_grid(lx, ly, dx, dy, Point3{Float64}, false);

subs = 7
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

function multi_solve(x)
  max_iter = 50
  atol = sqrt(eps(Float64))
  rtol = sqrt(eps(Float64))

  y = zeros(nv(sd))
  ares(y) = norm(Δ0 * y - x)
  rNorm = ares(y)
  ϵ = atol + rtol * rNorm
  count = 1
  while !(ares(y) ≤ ϵ || count ≥ max_iter)
    y = _multigrid_μ_cycle(y,x,md,gauss_seidel!,2)
    count += 1
  end
  y, (niter = count, )
end

mat = p2_d2_interpolation(sd)

left_wall_idxs = findall(p -> p[1] < dx, s[:point]);
right_wall_idxs = findall(p -> p[1] > lx - dx, s[:point]);

# For no change conditions in top and bottom cooling/heating elements
bottom_wall_idxs= findall(p -> p[2] == 0, s[:point]);
top_wall_idxs = findall(p -> p[2] == ly, s[:point]);

apply_tb_bc(x) = begin x[bottom_wall_idxs] .= 0; x[top_wall_idxs] .= 0; return x; end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => x -> begin
      y, _ = multi_solve(x)
      y .-= minimum(y)
    end
    :tb_bc => apply_tb_bc
    :interpolate => x -> mat * x
    _ => error("No operator $my_symbol found.")
  end
  return op
end

sim = eval(gensim(Porous_Convection))
f = sim(sd, generate, DiagonalHodge())

ΔT = 200.0

# Initial heat distribution
# Gaussian heat disturbance along with top and bottom elements
T_dist = MvNormal([lx/2.0, ly/2.0], [1/sqrt(2), 1/sqrt(2)])
T = [2 * ΔT * pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]
T[top_wall_idxs] .= -ΔT/2
T[bottom_wall_idxs] .= ΔT/2

# Measure the force of gravity in the downwards direction
accl_g = 9.81
# grav = SVector{3}([0.0, -accl_g, 0.0])
grav = SVector{3}([10.0, 20.0, 0.0])
grav = -accl_g * normalize(grav)
g = eval_constant_primal_form(sd, grav)
u₀ = ComponentArray(T=T, g=g)

# Physical constants
Ra = 1500
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

solvername = "MG"
filename = "Porous_Convection_$(solvername)_subs=$(subs)_dx=dy=$(dx)_Ra=$(Ra)"
save_dynamics("$(filename).mp4", length(soln.t))

iterations = Int64[]
for t in soln.t
  rho = codif_1 * wdg10(g, constants.αρ₀ * soln(t).T)
  _, stats = multi_solve(rho)
  push!(iterations, stats.niter)
end

save("$(filename)_histo.png", hist(iterations))
save("$(filename)_dots.png", plot(soln.t, iterations))

df = DataFrame(soln);
CSV.write("$(filename).csv", df)

# Can use to read out rows
# out_df = CSV.read("$(filename).csv", DataFrame)
# outT = Array(out_df[2, 2:end])
