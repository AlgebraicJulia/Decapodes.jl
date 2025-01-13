
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
import CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form

# Reference equations are at the link below, dual-time method is replaced with a Poisson equation solve
# https://pde-on-gpu.vaw.ethz.ch/lecture4/#solving_thermal_porous_convection_using_the_pseudo-transient_method

Porous_Convection = @decapode begin
  (λ_ρ₀Cp, αρ₀, k_ηf, ϕ)::Constant
  (P, T, Adv, bound_T, bound_Ṫ)::Form0
  (g, qD)::Form1

  bound_T == adiabatic(T)

  ρ == g ∧ (αρ₀ * bound_T)

  P == Δ⁻¹(δ(ρ))

  qD == -(k_ηf * (d(P) - ρ))

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

# fig = Figure();
# ax = CairoMakie.Axis(fig[1,1])
# msh = mesh!(ax, s)
# fig

Δ0 = Δ(0,sd);
fΔ0 = LinearAlgebra.factorize(Δ0);

# This matrix is for interpolation of PrimalForm2 to DualForm2, meant for use with primal-primal Lie
# TODO: To test, if there is a constant Form2, for example with just 1s, over primal triangles, primal vertex will have 1/(tri_area)
# This assume triangle area and orientations are the same. So essentially, inv_hdq_0 * mat * ones() = 1 ./ sd[:area]
mat = spzeros(nv(sd), ntriangles(sd))
for tri_id in triangles(sd)
  tri_area = sd[tri_id, :area]
  tri_orient = sd[tri_id, :tri_orientation]
  for dual_tri_id in tri_id:ntriangles(sd):nparts(sd, :DualTri)
    dual_tri_area = sd[dual_tri_id, :dual_area]
    dual_tri_orient = sd[dual_tri_id, :D_tri_orientation]

    weight = tri_orient * (dual_tri_area / tri_area)

    v = sd[sd[dual_tri_id, :D_∂e1], :D_∂v1]

    mat[v, tri_id] += weight

  end
end

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
      y .- minimum(y)
    end
    :adiabatic => x -> begin
      x[left_wall_idxs] .= x[next_left_wall_idxs]
      x[right_wall_idxs] .= x[next_right_wall_idxs]
      return x
    end
    :tb_bc => apply_tb_bc
    :interpolate => x -> mat * x
    _ => default_dec_matrix_generate(sd, my_symbol, hodge)
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

Ra = 1000
k_ηf = 1.0
αρ₀ = 1.0/accl_g
ϕ = 0.1
λ_ρ₀Cp = 1/Ra*(3.92*αρ₀*k_ηf*ΔT*ly/ϕ)
constants = (k_ηf = k_ηf, αρ₀ = αρ₀, ϕ = ϕ, λ_ρ₀Cp = λ_ρ₀Cp)

# Smallest time step in original simulation was 0.00019, largest around 0.00100, around 0.00050 from original implementation
# Only ran for 500 time steps, but adaptive time stepping means physical time simulated could vary
tₑ = 2

prob = ODEProblem(f, u₀, (0, tₑ), constants)
soln = solve(prob, Tsit5())

# For plotting Temperature, Pressure, Advection, and Diffusion

wdg10 = dec_wedge_product(Tuple{1, 0}, sd)
codif_1 = δ(1, sd)
d0 = dec_differential(0, sd)

hdg_1 = dec_hodge_star(1, sd)
dp_wdg_11 = dec_wedge_product_dp(Tuple{1,1}, sd)
inv_hdg_0 = dec_inv_hodge_star(0, sd)

function calculate_pressure(T, constants)
  boussinesq = wdg10(g, constants.αρ₀ * T)
  P = fΔ0 \ (codif_1 * boussinesq) # Pressure
end

function calculate_advection(T, P, constants)
  boussinesq = wdg10(g, constants.αρ₀ * T)

  darcy_flux = -(constants.k_ηf * (d0 * P - boussinesq))

  Adv = apply_tb_bc(-1/constants.ϕ *  inv_hdg_0 * mat * dp_wdg_11(hdg_1 * d0 * T, darcy_flux))
end

function calculate_diffusion(T, constants)
  Diff = apply_tb_bc(constants.λ_ρ₀Cp * Δ0 * T)
end

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)
  P = @lift(calculate_pressure($T, constants))
  Adv = @lift(calculate_advection($T, $P, constants))
  Diff = @lift(calculate_diffusion($T, constants))

  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Temperature at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s, color=T, colormap=:jet)
  Colorbar(f[1,2], msh_T)

  ax_P = CairoMakie.Axis(f[2,1], title = @lift("Pressure at Time $(round($time, digits=3))"))
  msh_P = mesh!(ax_P, s, color=P, colormap=:jet)
  Colorbar(f[2,2], msh_P)

  ax_Adv = CairoMakie.Axis(f[1,3], title = @lift("Advection at Time $(round($time, digits=3))"))
  msh_Adv = mesh!(ax_Adv, s, color=Adv, colormap=:jet)
  Colorbar(f[1,4], msh_Adv)

  ax_Diff = CairoMakie.Axis(f[2,3], title = @lift("Diffusion at Time $(round($time, digits=3))"))
  msh_Diff = mesh!(ax_Diff, s, color=Diff, colormap=:jet)
  Colorbar(f[2,4], msh_Diff)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("Porous_Convection.mp4", 120)
