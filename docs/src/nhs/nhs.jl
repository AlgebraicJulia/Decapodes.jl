# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics: Point3
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
Point3D = Point3{Float64};

@info "Creating Decapode"
momentum =  @decapode begin
  (v,V)::DualForm1
  f::Form0
  uˢ::DualForm1
  ∂tuˢ::DualForm1
  p::DualForm0
  b::DualForm0
  ĝ::DualForm1
  Fᵥ::DualForm1
  StressDivergence::DualForm1

  ∂ₜ(v) ==
    -ℒ₁(v,v) + 0.5*d(ι₁₁(v,v)) -
     d(ι₁₁(v,V)) + ι₁₂(v,d(V)) + ι₁₂(V,d(v)) -
     (f - ∘(d,⋆)(uˢ)) ∧ᵖᵈ₀₁ v -
     d(p) +
     b ∧ᵈᵈ₀₁ ĝ -
     StressDivergence +
     ∂tuˢ +
     Fᵥ
end

tracer_conservation = @decapode begin
  (c,C,F,FluxDivergence)::DualForm0
  (v,V)::DualForm1

  ∂ₜ(c) ==
    -1*ι₁₁(v,d(c)) -
    ι₁₁(V,d(c)) -
    ι₁₁(v,d(C)) -
    FluxDivergence +
    F
end

equation_of_state = @decapode begin
  (b,T,S)::DualForm0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end

isotropic_diffusivity = @decapode begin
  v::DualForm1
  c::DualForm0
  StressDivergence::DualForm1
  FluxDivergence::DualForm0
  (κ,nu)::Constant

  StressDivergence == nu*Δᵈ₁(v)
  FluxDivergence == κ*Δᵈ₀(c)
end

tracer_composition = @relation () begin
  # "The turbulence closure selected by the user determines the form of ... diffusive flux divergence"
  turbulence(FD,v,c)

  continuity(FD,v,c)
end

isotropic_tracer = apex(oapply(tracer_composition, [
  Open(isotropic_diffusivity, [:FluxDivergence, :v, :c]),
  Open(tracer_conservation,   [:FluxDivergence, :v, :c])]))

  nonhydrostatic_composition = @relation () begin
  # "The turbulence closure selected by the user determines the form of stress divergence"
  #   => Note that the StressDivergence term, SD, is shared by momentum and all the tracers.
  momentum(V, v, b, SD)

  # "Both T and S obey the tracer conservation equation"
  #   => Temperature and Salinity both receive a copy of the tracer physics.
  temperature(V, v, T, SD, nu)
  salinity(V, v, S, SD, nu)

  # "Buoyancy is determined from a linear equation of state"
  #   => The b term in momentum is that described by the equation of state here.
  eos(b, T, S)
end

isotropic_nonhydrostatic_buoyancy = apex(oapply(nonhydrostatic_composition, [
  Open(momentum,          [:V, :v, :b, :StressDivergence]),
  Open(isotropic_tracer,  [:continuity_V, :v, :c, :turbulence_StressDivergence, :turbulence_nu]),
  Open(isotropic_tracer,  [:continuity_V, :v, :c, :turbulence_StressDivergence, :turbulence_nu]),
  Open(equation_of_state, [:b, :T, :S])]));

@info "Creating Mesh"
# This is a torus with resolution of its dual mesh similar to that
# used by Oceananigans (explicitly represented as a torus, not as a
# square with periodic boundary conditions!)
download("https://cise.ufl.edu/~luke.morris/torus.obj", "torus.obj")
s = EmbeddedDeltaSet2D("torus.obj")
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Barycenter())

@info "Generating Simulation"
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    _ => default_dec_matrix_generate(sd, my_symbol, hodge)
  end
  return op
end

open("nhs_sim.jl", "w") do f
  write(f, string(gensim(isotropic_nonhydrostatic_buoyancy)))
end
sim = include("nhs_sim.jl")

@info "Creating Initial Conditions"
fₘ = sim(sd, generate)

S = map(sd[sd[:tri_center], :dual_point]) do (_,_,_)
  0.0
end
T = map(sd[sd[:tri_center], :dual_point]) do (_,_,_)
  0.0
end
p = map(sd[sd[:tri_center], :dual_point]) do (_,_,_)
  0.0
end
f = zeros(nv(sd))
Fₛ = zeros(ntriangles(sd))
Fₜ = zeros(ntriangles(sd))
Cₛ = zeros(ntriangles(sd))
Cₜ = zeros(ntriangles(sd))
V = zeros(ne(sd))
v = rand(ne(sd)) * 1e-8
ĝ = ♭(sd, DualVectorField(fill(Point3D(0,0,0), ntriangles(sd)))).data
Fᵥ = zeros(ne(sd))
qₛ = zeros(ne(sd))
qₜ = zeros(ne(sd))
uˢ = zeros(ne(sd))
∂tuˢ = zeros(ne(sd))

u₀ = ComponentArray(
  v = v,
  V = V,
  momentum_f = f,
  momentum_uˢ = uˢ,
  momentum_∂tuˢ = ∂tuˢ,
  momentum_p = p,
  momentum_ĝ = ĝ,
  momentum_Fᵥ = Fᵥ,
  T = T,
  temperature_continuity_C = Cₜ,
  temperature_continuity_F = Fₜ,
  S = S,
  salinity_continuity_C = Cₛ,
  salinity_continuity_F = Fₛ)

gᶜ = 9.81
α = 2e-3
β = 5e-4
constants_and_parameters = (
  temperature_turbulence_κ = 0.0,
  salinity_turbulence_κ = 0.0,
  nu = 1e-5,
  eos_g = gᶜ,
  eos_α = α,
  eos_β = β);

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Vern7())
soln.retcode != :Unstable || error("Solver was not stable")

tₑ = 50

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Vern7(), force_dtmin=true, dtmax=0.2)
@show soln.retcode
@info("Done")

ihs0 = dec_inv_hodge_star(Val{0}, sd, GeometricHodge())
dd1 = dec_dual_derivative(1, sd)
♯_m = ♯_mat(sd, LLSDDSharp())
using LinearAlgebra: norm
function vorticity(α)
  ihs0*dd1*α
end
function speed(α)
  norm.(♯_m * α)
end

function save_vorticity(is_2d=false)
  frames = 200
  time = Observable(0.0)
  fig = Figure(title = @lift("Vorticity at $($time)"))
  ax = is_2d ?
    CairoMakie.Axis(fig[1,1]) :
    LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.mesh!(ax, s,
    color=@lift(vorticity(soln($time).v)),
    colorrange=extrema(vorticity(soln(tₑ).v)).*.9,
    colormap=:jet)

  Colorbar(fig[1,2], msh)
  record(fig, "vorticity.gif", range(0.0, tₑ; length=frames); framerate = 20) do t
    time[] = t
  end
end
save_vorticity(false)

function save_speed(is_2d=false) frames = 200
  time = Observable(0.0)
  fig = Figure(title = @lift("Speed at $($time)"))
  ax = is_2d ?
    CairoMakie.Axis(fig[1,1]) :
    LScene(fig[1,1], scenekw=(lights=[],))
  msh = CairoMakie.scatter!(ax, sd[sd[:tri_center], :dual_point],
    color=@lift(speed(soln($time).v)),
    colorrange=extrema(speed(soln(tₑ).v)).*.9,
    colormap=:jet,
    markersize=5)

  Colorbar(fig[1,2], msh)
  record(fig, "speed.gif", range(0.0, tₑ; length=frames); framerate = 20) do t
    time[] = t
  end
end
save_speed(false)
