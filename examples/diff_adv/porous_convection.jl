
using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using Catlab
using LinearAlgebra
using CairoMakie
using GeometryBasics: Point2, Point3
using Distributions
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using SparseArrays
using StaticArrays
using StatsBase
import CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form

# Reference equations are at the link below, dual-time method is replaced with a Poisson equation solve
# https://pde-on-gpu.vaw.ethz.ch/lecture4/#solving_thermal_porous_convection_using_the_pseudo-transient_method

Porous_Convection = @decapode begin
  (λ_ρ₀Cp, αρ₀, k_ηf, ϕ)::Constant
  (P, T, Adv, bound_T, bound_Ṫ)::Form0
  (g, qD)::Form1

  bound_T == adiabatic(T)
  # Darcy flux
  ρ == g ∧ (αρ₀ * bound_T)
  P == Δ⁻¹(δ(ρ))
  qD == -k_ηf * (d(P) - ρ)

  Adv == ⋆(interpolate(∧ᵈᵖ₁₁(⋆(d(bound_T)), qD)))
  Ṫ == -1/ϕ * Adv + λ_ρ₀Cp * Δ(bound_T)

  bound_Ṫ == tb_bc(Ṫ)

  ∂ₜ(T) == bound_Ṫ
end
infer_types!(Porous_Convection)
resolve_overloads!(Porous_Convection)
# to_graphviz(Porous_Convection)

lx, ly = 40.0, 20.0
dx = dy = 0.4
s = triangulated_grid(lx, ly, dx, dy, Point3{Float64});
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point2{Float64}}(s);
subdivide_duals!(sd, Circumcenter());

Δ0 = Δ(0,sd);
fΔ0 = LinearAlgebra.factorize(Δ0);

mat = p2_d2_interpolation(sd)

left_wall_idxs = findall(p -> p[1] < dx, s[:point]);
right_wall_idxs = findall(p -> p[1] > lx - dx, s[:point]);

# For adiabatic conditions: https://en.wikipedia.org/wiki/Adiabatic_process
next_left_wall_idxs = left_wall_idxs .+ 1;
next_right_wall_idxs = right_wall_idxs .- 1;

# For no change conditions in top and bottom cooling/heating elements
bottom_wall_idxs= findall(p -> p[2] == 0, s[:point]);
top_wall_idxs = findall(p -> p[2] == ly, s[:point]);

apply_tb_bc(x) = begin x[bottom_wall_idxs] .= 0; x[top_wall_idxs] .= 0; return x; end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => x -> begin
      y = fΔ0 \ x
      y .-= minimum(y)
    end
    :adiabatic => x -> begin
      x[left_wall_idxs] .= x[next_left_wall_idxs]
      x[right_wall_idxs] .= x[next_right_wall_idxs]
      return x
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

T_dist = MvNormal([lx/2.0, ly/2.0], [1/sqrt(2), 1/sqrt(2)])
T = [2 * ΔT * pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]
T[top_wall_idxs] .= -ΔT/2
T[bottom_wall_idxs] .= ΔT/2

# Measure the force of gravity in the downwards direction
accl_g = 9.81
grav = SVector{3}([0.0, -accl_g, 0.0])
g = eval_constant_primal_form(sd, grav)
u₀ = ComponentArray(T=T, g=g)

Ra = 750
k_ηf = 1.0
αρ₀ = (1.0/accl_g)
ϕ = 0.1
λ_ρ₀Cp = 1/Ra*(accl_g*αρ₀*k_ηf*ΔT*ly/ϕ)
constants = (k_ηf = k_ηf, αρ₀ = αρ₀, ϕ = ϕ, λ_ρ₀Cp = λ_ρ₀Cp)

# Smallest time step in original simulation was 0.00019, largest around 0.00100, around 0.00050 from original implementation
# Only ran for 500 time steps, but adaptive time stepping means physical time simulated could vary
tₑ = 0.7

prob = ODEProblem(f, u₀, (0, tₑ), constants)
soln = solve(prob, Tsit5(); saveat = 0.005)

# For plotting Temperature, Pressure, Advection, and Diffusion

wdg10 = dec_wedge_product(Tuple{1, 0}, sd)
codif_1 = δ(1, sd)
d0 = dec_differential(0, sd)

hdg_1 = dec_hodge_star(1, sd)
dp_wdg_11 = dec_wedge_product_dp(Tuple{1,1}, sd)
inv_hdg_0 = dec_inv_hodge_star(0, sd)

function calculate_pressure(T, constants)
  fΔ0 \ (codif_1 * wdg10(g, constants.αρ₀ * T))
end

function calculate_advection(T, P, constants)
  darcy_flux = -constants.k_ηf * (d0 * P - wdg10(g, constants.αρ₀ * T))
  apply_tb_bc(-1/constants.ϕ *  inv_hdg_0 * mat * dp_wdg_11(hdg_1 * d0 * T, darcy_flux))
end

function calculate_diffusion(T, constants)
  apply_tb_bc(constants.λ_ρ₀Cp * Δ0 * T)
end

function compute_colorranges(length)
  values = hcat(map(range(0, soln.t[end], length=length)) do t
    T = soln(t).T
    P = calculate_pressure(T, constants)
    Adv = calculate_advection(T, P, constants)
    Diff = calculate_diffusion(T, constants)

    [minimum(P), maximum(P), minimum(Adv), maximum(Adv), minimum(Diff), maximum(Diff)]
  end...)

  # Minimum, Maximum
  percentile(values[1, :], 90), percentile(values[2, :], 90), # Pressure
  percentile(values[3, :], 90), percentile(values[4, :], 90), # Advection
  percentile(values[5, :], 75), percentile(values[6, :], 75) # Diffusion
end

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)
  P = @lift(calculate_pressure($T, constants))
  Adv = @lift(calculate_advection($T, $P, constants))
  Diff = @lift(calculate_diffusion($T, constants))

  colorranges = compute_colorranges(video_length)
  P_range = colorranges[1], colorranges[2]
  Adv_range = colorranges[3], colorranges[4]
  Diff_range = colorranges[5], colorranges[6]

  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Temperature at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s; color=T, colormap=:jet, colorrange=(-ΔT/2, ΔT/2))
  Colorbar(f[1,2], msh_T)

  ax_P = CairoMakie.Axis(f[2,1], title = @lift("Pressure at Time $(round($time, digits=3))"))
  msh_P = mesh!(ax_P, s; color=P, colormap=:jet, colorrange=P_range)
  Colorbar(f[2,2], msh_P)

  ax_Adv = CairoMakie.Axis(f[1,3], title = @lift("Advection at Time $(round($time, digits=3))"))
  msh_Adv = mesh!(ax_Adv, s; color=Adv, colormap=:jet, colorrange=Adv_range)
  Colorbar(f[1,4], msh_Adv)

  ax_Diff = CairoMakie.Axis(f[2,3], title = @lift("Diffusion at Time $(round($time, digits=3))"))
  msh_Diff = mesh!(ax_Diff, s; color=Diff, colormap=:jet, colorrange=Diff_range)
  Colorbar(f[2,4], msh_Diff)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("Porous_Convection_New.mp4", 120)
