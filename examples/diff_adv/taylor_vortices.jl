using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using CairoMakie
using ComponentArrays
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using CoordRefSystems
using SparseArrays

using LoggingExtras
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

euler_equations = @decapode begin
  du::DualForm2
  u::DualForm1
  ψ::Form0
  ω::Form0
  μ::Constant

  ψ == Δ⁻¹(⋆(du))
  u == ⋆(d(ψ))
  v == ♭♯(u)

  ω == ⋆(d(u) + dᵦ(v))

  ∂ₜ(du) ==  -(∘(⋆₁, d̃₁)(∧₁₀(v, ω))) + μ*d(⋆(d(ω)))
end

const RADIUS = 1.0
s = loadmesh(Icosphere(7, RADIUS));
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter());

Δ0 = Δ(0,sd);
fΔ0 = factorize(Δ0);
dual_d0 = dec_dual_derivative(0,sd);
dual_d1 = dec_dual_derivative(1,sd);
dᵦ = 0.5 * abs.(dual_d1) * spdiagm(dual_d0 * ones(ntriangles(sd)));
function generate(s, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => x -> begin
      y = fΔ0 \ x
      y .- minimum(y)
    end
    :dᵦ => x -> dᵦ * x
  _ => error("Unmatched operator $my_symbol")
  end
  return op
end;

sim = evalsim(euler_equations)
fₘ = sim(sd, generate, DiagonalHodge());

struct TaylorVortexParams
  G::Real
  a::Real
end

function great_circle_dist(pnt1::Point3D, pnt2::Point3D)
  RADIUS * acos(dot(pnt1,pnt2))
end

function taylor_vortex(pnt::Point3D, cntr::Point3D, p::TaylorVortexParams)
  gcd = great_circle_dist(pnt,cntr)
  (p.G/p.a) * (2 - (gcd/p.a)^2) * exp(0.5 * (1 - (gcd/p.a)^2))
end

taylor_vortex(sd::HasDeltaSet, cntr::Point3D, p::TaylorVortexParams) =
  map(x -> taylor_vortex(x, cntr, p), point(sd))

function vort_ring(lat, n_vorts, p::TaylorVortexParams, formula)
  sum(map(x -> formula(sd, x, p), ring_centers(lat, n_vorts)))
end

function ring_centers(lat, n)
  ϕs = range(0.0, 2π; length=n+1)[1:n]
  map(ϕs) do ϕ
    v_sph = Spherical(RADIUS, lat, ϕ)
    v_crt = convert(Cartesian, v_sph)
    Point3D(v_crt.x.val, v_crt.y.val, v_crt.z.val)
  end
end

X = vort_ring(0.2, 2, TaylorVortexParams(0.5, 0.1), taylor_vortex)

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=X, colormap=Reverse(:redsblues))
Colorbar(fig[1,2], msh)
fig

tₑ = 10.0

hdg_0 = dec_hodge_star(0,sd,DiagonalHodge())
u₀ = ComponentArray(du = hdg_0*X)

prob = ODEProblem(fₘ, u₀, (0, tₑ), (μ=1e-3,))
soln = solve(prob, Tsit5(), dtmax = 0.01, saveat=tₑ/10.0, dense=false, progress=true, progress_steps=1);

inv_hdg_0 = dec_inv_hodge_star(0,sd,DiagonalHodge())
fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=inv_hdg_0 * soln.u[7].du, colormap=Reverse(:redsblues))
Colorbar(fig[1,2], msh)
fig
