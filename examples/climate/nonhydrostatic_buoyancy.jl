# This is a 50 line implementation of Oceananigans.jl's NonhydrostaticModel.

#######################
# Import Dependencies #
#######################

# AlgebraicJulia Dependencies
using Catlab
using CombinatorialSpaces
using Decapodes

# External Dependencies
using GeometryBasics: Point3
using GLMakie
using JLD2
using LinearAlgebra
using MultiScaleArrays
using OrdinaryDiffEq
Point3D = Point3{Float64}

####################
# Define the model #
####################

momentum = @decapode begin
  (f,b)::Form0
  (v,V,g,Fᵥ,uˢ)::Form1
  τ::Form2
  v̇ == ∂ₜ(v)

  v̇ == -1 * L(v,v) - L(V,v) - L(v,V) -
       f∧v - ∘(⋆,d,⋆)(uˢ)∧v - d(p) + b∧g - ∘(⋆,d,⋆)(τ) + dtuˢ + Fᵥ
end
momentum = expand_operators(momentum)
to_graphviz(momentum)

tracer = @decapode begin
  (c,F)::Form0
  (v,V,q)::Form1
  ċ == ∂ₜ(c)

  ċ == -1*⋆(L(v,⋆(c))) - ⋆(L(V,⋆(c))) - ⋆(L(v,⋆(c))) - ∘(⋆,d,⋆)(q) + F
end
to_graphviz(tracer)

equation_of_state = @decapode begin
  (b,T,S)::Form0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end
to_graphviz(equation_of_state)

buoyancy_composition_diagram = @relation () begin
  momentum(V, v, b)

  # "Both T and S obey the tracer conservation equation"
  temperature(V, v, T)
  salinity(V, v, S)

  # "Buoyancy is determined from a linear equation of state"
  eos(b, T, S)
end
to_graphviz(buoyancy_composition_diagram, box_labels=:name, junction_labels=:variable, prog="fdp", graph_attrs=Dict(["sep" => "1.5"]))

buoyancy_cospan = oapply(buoyancy_composition_diagram, [
  Open(momentum,          [:V, :v, :b]),
  Open(tracer,            [:V, :v, :c]),
  Open(tracer,            [:V, :v, :c]),
  Open(equation_of_state, [:b, :T, :S])])

buoyancy = apex(buoyancy_cospan)
to_graphviz(buoyancy)

buoyancy = expand_operators(buoyancy)
infer_types!(buoyancy)
resolve_overloads!(buoyancy)
to_graphviz(buoyancy)

s′ = loadmesh(Rectangle_30x10())
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())

include("dec_operators.jl")
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :L₁ => (x,y) -> L₁′(x,y,sd,hodge)
    :L₂ᵈ => (x,y) -> lie_derivative_flat(2, sd, x, y)
    _ => default_dec_generate(sd, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end

sim = eval(gensim(buoyancy))
fₘ = sim(s, generate)

S = zeros(nv(s))
T = zeros(nv(s))
p = zeros(nv(s))
f = zeros(nv(s))
Fₛ = zeros(nv(s))
Fₜ = zeros(nv(s))
f = zeros(nv(s))

V = zeros(ne(s))
v = zeros(ne(s))
g = ♭(sd, DualVectorField(fill(Point3D(0,1,0), ntriangles(s)))).data
Fᵥ = zeros(ne(s))
qₛ = zeros(ne(s))
qₜ = zeros(ne(s))
uˢ = zeros(ne(s))
dtuˢ = zeros(ne(s))

τ = zeros(ntriangles(s))

u₀ = construct(PhysicsState, [
  VectorForm(f),
  VectorForm(v),
  VectorForm(V),
  VectorForm(g),
  VectorForm(Fᵥ),
  VectorForm(uˢ),
  VectorForm(τ),
  VectorForm(dtuˢ),
  VectorForm(p),
  VectorForm(T),
  VectorForm(Fₜ),
  VectorForm(qₜ),
  VectorForm(S),
  VectorForm(Fₛ),
  VectorForm(qₛ)], Float64[], [
  :momentum_f,
  :v,
  :V,
  :momentum_g,
  :momentum_Fᵥ,
  :momentum_uˢ,
  :momentum_τ,
  :momentum_dtuˢ,
  :momentum_p,
  :T,
  :temperature_F,
  :temperature_q,
  :S,
  :salinity_F,
  :salinity_q])

gᶜ = 0
α = 0
β = 0
constants_and_parameters = (
  eos_g = gᶜ,
  eos_α = α,
  eos_β = β)

tₑ = 1e9

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

begin end

## Track a single tracer.
#single_tracer_composition_diagram = @relation () begin
#  momentum(V, v)
#  tracer(V, v)
#end
#to_graphviz(single_tracer_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")
#
#single_tracer_cospan = oapply(single_tracer_composition_diagram,
#  [Open(momentum, [:V, :v]),
#  Open(tracer, [:V, :v])])
#
#single_tracer = apex(single_tracer_cospan)
#to_graphviz(single_tracer)
#
#single_tracer = expand_operators(single_tracer)
#infer_types!(single_tracer)
#resolve_overloads!(single_tracer)
#to_graphviz(single_tracer)
#
#
## Track multiple tracers.
#triple_tracer_composition_diagram = @relation () begin
#  momentum(V, v)
#
#  mercury(V, v)
#  phosphate(V, v)
#  oil(V, v)
#end
#to_graphviz(triple_tracer_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")
#
#triple_tracer_cospan = oapply(triple_tracer_composition_diagram,
#  [Open(momentum, [:V, :v]),
#  Open(tracer, [:V, :v]),
#  Open(tracer, [:V, :v]),
#  Open(tracer, [:V, :v])])
#
#triple_tracer = apex(triple_tracer_cospan)
#to_graphviz(triple_tracer)
#
#triple_tracer = expand_operators(triple_tracer)
#infer_types!(triple_tracer)
#resolve_overloads!(triple_tracer)
#to_graphviz(triple_tracer)
#
#
## The buoyancy model uses b as a tracer.
#buoyancy_composition_diagram = @relation () begin
#  momentum(V, v, b)
#
#  # "b obeys the tracer equation"
#  tracer(V, v, b)
#end
#to_graphviz(buoyancy_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")
#
#buoyancy_cospan = oapply(buoyancy_composition_diagram,
#  [Open(momentum, [:V, :v, :b]),
#  Open(tracer, [:V, :v, :c])])
#
#buoyancy = apex(buoyancy_cospan)
#to_graphviz(buoyancy)
#
#buoyancy = expand_operators(buoyancy)
#infer_types!(buoyancy)
#resolve_overloads!(buoyancy)
#to_graphviz(buoyancy)
