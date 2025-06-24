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
# incompressible_ns = @decapode begin
#   du::DualForm2
#   u::DualForm1
#   ψ::Form0
#   ω::Form0
#   μ::Constant

#   ψ == Δ⁻¹(⋆(du))
#   u == ⋆(d(ψ))
#   v == ♭♯(u)

#   ω == ⋆(d(u) + dᵦ(v))

#   ∂ₜ(du) == -(d(⋆(∧(v, ω)))) + μ*d(⋆(d(ω)))
# end

incompressible_ns = @decapode begin
  ψ::Form0
  u::DualForm1
  v::Form1
  ω::Form0
  μ::Constant

  u == ⋆(d(ψ))
  v == ♭♯(u)

  ω == ⋆(d(u) + dᵦ(v))

  ∂ₜ(ψ) == dsdinv(-(d(⋆(∧(v, ω)))) + μ*d(⋆(d(ω))))
end

s = triangulated_grid(2π,2π,0.02,0.02,Point3d,false)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter());

dual_d0 = dec_dual_derivative(0,sd);
dual_d1 = dec_dual_derivative(1,sd);
dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));
hdg_1 = dec_hodge_star(1, sd, DiagonalHodge())
d0 = dec_differential(0, sd)
dsd = dual_d1 * hdg_1 * d0
dsdinv = factorize(dsd)

function generate(s, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => x -> begin
      y = fΔ0 \ x
      y .- minimum(y)
    end
    :dsdinv => x -> dsdinv \ x
    :dᵦ => x -> dᵦ * x
  _ => error("Unmatched operator $my_symbol")
  end
  return op
end;

sim = evalsim(incompressible_ns)
fₘ = sim(sd, generate, DiagonalHodge());

struct TaylorVortexParams
  G::Real
  a::Real
end

function taylor_vortex(pnt::Point3d, cntr::Point3d, p::TaylorVortexParams)
  r = norm(pnt .- cntr)
  (p.G/p.a) * (2 - (r/p.a)^2) * exp(0.5 * (1 - (r/p.a)^2))
end

taylor_vortex(sd::HasDeltaSet, cntr::Point3d, p::TaylorVortexParams) =
  map(x -> taylor_vortex(x, cntr, p), point(sd))

function vort_ring(p::TaylorVortexParams, formula)
  sum(map(x -> formula(sd, x, p), [Point3d(π-0.4, π, 0.0), Point3d(π+0.4, π, 0.0)]))
end

ω = vort_ring(TaylorVortexParams(1.0, 0.3), taylor_vortex)

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=ω, colormap=:jet)
Colorbar(fig[1,2], msh)
display(fig)

tₑ = 10.0

Δ0 = Δ(0,sd);
fΔ0 = factorize(Δ0);

ψ = fΔ0 \ ω
ψ .= ψ .- minimum(ψ)

hdg_0 = dec_hodge_star(0,sd,DiagonalHodge())
u₀ = ComponentArray(ψ = ψ)

μ = 0
prob = ODEProblem(fₘ, u₀, (0, tₑ), (μ=μ,))
soln = solve(prob, Tsit5(), dtmax = 0.01, saveat=tₑ/50.0, progress=true, progress_steps=1);

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=Δ0*soln.u[end].ψ, colormap=:jet)
Colorbar(fig[1,2], msh)
display(fig)

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  w = @lift(Δ0 * soln($time).ψ)

  f = Figure()

  ax_w = CairoMakie.Axis(f[1,1], title = @lift("Vorticity at Time $(round($time, digits=2)) with μ=$(μ)"))
  msh_w = mesh!(ax_w, s; color=w, colormap=:jet, colorrange=extrema(ω))
  Colorbar(f[1,2], msh_w)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("Taylor_Vortices_streamfunction_mu=$(μ).mp4", 30)
