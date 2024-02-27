# This is a discrete exterior calculus implementation of Oceananigans.jl's `NonhydrostaticModel`.

#######################
# Import Dependencies #
#######################

# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using CombinatorialSpaces: interior_product_dd, ℒ_dd
using Decapodes
using DiagrammaticEquations

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics: Point3
using JLD2
using LinearAlgebra
using Logging: global_logger
using MLStyle
using OrdinaryDiffEq
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
Point3D = Point3{Float64}

####################
# Define the model #
####################

# Equation 1: "The momentum conservation equation" from https://clima.github.io/OceananigansDocumentation/stable/physics/nonhydrostatic_model/#The-momentum-conservation-equation
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
    -ℒ(v,v) + 0.5*d(ι₁₁(v,v)) -
     d(ι₁₁(v,V)) + ι₁₂(v,d(V)) + ι₁₂(V,d(v)) -
     (f - ∘(d,⋆)(uˢ)) ∧ᵖ v -
     d(p) +
     b ∧ᵈ ĝ -
     StressDivergence +
     ∂tuˢ +
     Fᵥ
end

# Why write "StressDivergence" instead of ∇⋅τ?
# According to this docs page:
# https://clima.github.io/OceananigansDocumentation/stable/physics/turbulence_closures/
# , the user simply selects what model to insert in place of the term ∇⋅τ.
# For example, in the isotropic case, Oceananigans.jl rewrites with:
# ∇⋅τ = nuΔv.
# Thus, we write StressDivergence, and replace this term with a choice of
# "turbulence closure" model. Using the "constant isotropic diffusivity" case,
# we can operate purely in terms of scalar-valued forms.

# Equation 2: "The tracer conservation equation" from https://clima.github.io/OceananigansDocumentation/stable/physics/nonhydrostatic_model/#The-tracer-conservation-equation
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

# Equation 2: "Linear equation of state" of seawater buoyancy from https://clima.github.io/OceananigansDocumentation/stable/physics/buoyancy_and_equations_of_state/#Linear-equation-of-state
equation_of_state = @decapode begin
  (b,T,S)::DualForm0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end

# Equation 2: "Constant isotropic diffusivity" from https://clima.github.io/OceananigansDocumentation/stable/physics/turbulence_closures/#Constant-isotropic-diffusivity
isotropic_diffusivity = @decapode begin
  v::DualForm1
  c::DualForm0
  StressDivergence::DualForm1
  FluxDivergence::DualForm0
  (κ,nu)::Constant

  StressDivergence == nu*Δᵈ₁(v)
  FluxDivergence == κ*Δᵈ₀(c)
end

##################
# Compose Models #
##################

"""
Decapodes composition is formally known as an "operad algebra". That means that we don't have to encode our composition in a single UWD and then apply it. Rather, we can define several UWDs, compose those, and then apply those. Of course, since the output of oapply is another Decapode, we could perform an intermediate oapply, if that is convenient.

Besides it being convenient to break apart large UWDs into component UWDs, this hierarchical composition can enforce rules on our physical quantities.

For example:
1. We want all the tracers (salinity, temperature, etc.) in our physics to obey the same conservation equation.
2. We want them to obey the same "turbulence closure", which affects their flux-divergence term.
3. At the same time, a choice of turbulence closure doesn't just affect (each of) the flux-divergence terms, it also constrains which stress-divergence is physically valid in the momentum equation.

We will use our operad algebra to guarantee model compatibility and physical consistency, guarantees that would be burdensome to fit into a one-off type system.
"""

# Specify the equations that a tracer obeys:
tracer_composition = @relation () begin
  # "The turbulence closure selected by the user determines the form of ... diffusive flux divergence"
  turbulence(FD,v,c)

  continuity(FD,v,c)
end

# Let's "lock in" isotropic diffusivity by doing an intermediate oapply.
isotropic_tracer = apex(oapply(tracer_composition, [
  Open(isotropic_diffusivity, [:FluxDivergence, :v, :c]),
  Open(tracer_conservation,   [:FluxDivergence, :v, :c])]))

# Use this building-block tracer physics at the next level:

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
  Open(equation_of_state, [:b, :T, :S])]))

#################
# Define a Mesh #
#################

# This is the sphere with the same surface area as the Oceananigans
# example torus:
#s = loadmesh(Icosphere(5, (3*π)^(1/3)))

# This is a torus with resolution of its dual mesh similar to that
# used by Oceananigans (explicitly represented as a torus, not as a
# square with periodic boundary conditions!)
#download("https://cise.ufl.edu/~luke.morris/torus.obj", "torus.obj")
s = EmbeddedDeltaSet2D("torus.obj")
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Barycenter())
xmax = maximum(x -> x[1], point(s))
zmax = maximum(x -> x[2], point(s))
#wireframe(s)

#################################
# Define Differential Operators #
#################################

# TODO: Provide these functions by default.
i11 = interior_product_dd(Tuple{1,1}, sd);
i12 = interior_product_dd(Tuple{1,2}, sd);
lie11 = ℒ_dd(Tuple{1,1}, sd);
Λᵖ = dec_wedge_product_pd(Tuple{0,1}, sd);
Λᵈ = dec_wedge_product_dd(Tuple{0,1}, sd);
# TODO: Upstream the dual 0 Laplacian.
dd0 = dec_dual_derivative(0, sd);
ihs1 = dec_inv_hodge_star(1, sd, GeometricHodge());
d1 = dec_differential(1,sd);
hs2 = dec_hodge_star(2, sd, GeometricHodge());
# TODO: Upstream the dual 1 Laplacian.
dd1 = dec_dual_derivative(1, sd);
ihs0 = dec_inv_hodge_star(0, sd, GeometricHodge());
d0 = dec_differential(0,sd);
hs1 = dec_hodge_star(1, sd, GeometricHodge());
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :ℒ => lie11
    :ι₁₁ => i11
    :ι₁₂ => i12
    :∧ᵖ => Λᵖ
    :∧ᵈ => Λᵈ
    :Δᵈ₀ => x -> hs2 * d1 * ihs1(dd0 * x)
    :Δᵈ₁ => x -> begin
      hs1 * d0 * ihs0 * dd1 * x +
      dd0 * hs2 * d1 * ihs1(x)
    end
    _ => default_dec_generate(sd, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end

#######################
# Generate Simulation #
#######################

open("nhs.jl", "w") do f
  write(f, string(gensim(expand_operators(isotropic_nonhydrostatic_buoyancy))))
end
sim = include("nhs.jl")
fₘ = sim(sd, generate)

#################################
# Define Differential Operators #
#################################

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
#v = zeros(ne(sd))
#v = rand(ne(sd)) *1e-5
v = rand(ne(sd)) *1e-8
#ĝ = ♭(sd, DualVectorField(fill(Point3D(0,1,0), ntriangles(sd)))).data
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
  nu = 0.0,
  eos_g = gᶜ,
  eos_α = α,
  eos_β = β)
  

tₑ = 5e4

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5(), dtmin=1e-3, force_dtmin=true)
soln.retcode != :Unstable || error("Solver was not stable")

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5(), dtmin=1e-16, force_dtmin=true, progress=true)
@show soln.retcode
@info("Done")

# Vorticity:
mesh(s,  
        color=ihs0*dd1*soln(tₑ).v,
        colormap=:jet)

# Speed:

# Create a gif
begin
  frames = 100
  #fig, ax, ob = GLMakie.mesh(s′, color=soln(0).T, colormap=:jet, colorrange=extrema(soln(tₑ).h))
  fig, ax, ob = GLMakie.mesh(s′, color=soln(0).T, colormap=:jet, colorrange=extrema(soln(1.5).T))
  Colorbar(fig[1,2], ob)
  #record(fig, "oceananigans.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
  record(fig, "oceananigans.gif", range(0.0, 1.5; length=frames); framerate = 30) do t
    ob.color = soln(t).T
  end
end

begin end

# Track a single tracer.
single_tracer_composition_diagram = @relation () begin
  momentum(V, v)
  tracer(V, v)
end
to_graphviz(single_tracer_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")

single_tracer_cospan = oapply(single_tracer_composition_diagram,
  [Open(momentum, [:V, :v]),
  Open(tracer, [:V, :v])])

single_tracer = apex(single_tracer_cospan)
to_graphviz(single_tracer)

single_tracer = expand_operators(single_tracer)
infer_types!(single_tracer)
resolve_overloads!(single_tracer)
to_graphviz(single_tracer)


# Track multiple tracers.
triple_tracer_composition_diagram = @relation () begin
  momentum(V, v)

  mercury(V, v)
  phosphate(V, v)
  oil(V, v)
end
to_graphviz(triple_tracer_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")

triple_tracer_cospan = oapply(triple_tracer_composition_diagram,
  [Open(momentum, [:V, :v]),
  Open(tracer, [:V, :v]),
  Open(tracer, [:V, :v]),
  Open(tracer, [:V, :v])])

triple_tracer = apex(triple_tracer_cospan)
to_graphviz(triple_tracer)

triple_tracer = expand_operators(triple_tracer)
infer_types!(triple_tracer)
resolve_overloads!(triple_tracer)
to_graphviz(triple_tracer)

