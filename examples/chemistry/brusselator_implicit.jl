using ACSets
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes
using MLStyle
using OrdinaryDiffEq
using Symbolics
using Krylov
using LinearSolve
using LinearAlgebra
using CairoMakie
using ComponentArrays
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

Brusselator = @decapode begin
  (U, V)::Form0
  U2V::Form0
  (U̇, V̇)::Form0

  (α)::Constant
  F::Parameter

  U2V == U .^ 2 .* V

  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(V))
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

x = 1
dx = 0.005
s = triangulated_grid(x,x,dx,dx,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());

U = map(sd[:point]) do (_,y)
  22 * (y * (1-y))^(3/2)
end

V = map(sd[:point]) do (x,_)
  27 * (x * (1-x))^(3/2)
end

F₁ = map(sd[:point]) do (x,y)
  (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end

F₂ = zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₁ : F₂)

sim = evalsim(Brusselator, can_prealloc=false)
fₘ = sim(sd, nothing, DiagonalHodge())

u₀ = ComponentArray(U=U, V=V)

du₀ = copy(u₀)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> fₘ(du, u, constants_and_parameters, 0.0), du₀, u₀)

f = ODEFunction(fₘ; jac_prototype = float.(jac_sparsity))
# tₑ = 11.5
tₑ = 20
@info("Solving")
prob = ODEProblem(f, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, FBDF(linsolve = KLUFactorization()), progress=true, progress_steps=1, saveat=0.1);
@info("Done")
soln.stats

function save_dynamics(save_file_name)
  time = Observable(0.0)
  u = @lift(soln($time).U)
  f = Figure()
  ax = CairoMakie.Axis(f[1,1], title = @lift("Brusselator U Concentration at Time $($time)"))
  gmsh = mesh!(ax, s, color=u, colormap=:jet, colorrange=extrema(soln(0.0).U))
  Colorbar(f[1,2], gmsh)
  timestamps = range(0, tₑ, step=1e-1)
  record(f, save_file_name, timestamps; framerate = 30) do t
    time[] = t
  end
end
save_dynamics("brusselator_U.gif")
