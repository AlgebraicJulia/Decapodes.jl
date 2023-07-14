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
  U::Parameter

  v̇ == ∂ₜ(v)
  uˢ̇ == ∂ₜ(uˢ)

  v̇ == -1 * L(v,v) - L(V,v) - L(v,V) -
       f∧v - ∘(⋆,d,⋆)(uˢ)∧v - d(p) + b∧g - ∘(⋆,d,⋆)(τ) + uˢ̇ + Fᵥ

  uˢ̇ == force(U)
end
momentum = expand_operators(momentum)
to_graphviz(momentum)

tracer = @decapode begin
  (c,C,F,c_up)::Form0
  (v,V,q)::Form1

  c_up == -1*⋆(L(v,⋆(c))) - ⋆(L(V,⋆(c))) - ⋆(L(v,⋆(C))) - ∘(⋆,d,⋆)(q) + F
end
to_graphviz(tracer)

equation_of_state = @decapode begin
  (b,T,S)::Form0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end
to_graphviz(equation_of_state)

boundary_conditions = @decapode begin
  (S,T)::Form0
  (Ṡ,T_up)::Form0
  Ṫ == ∂ₜ(T)
  Ṡ == ∂ₜ(S)

  Ṫ == ∂_spatial(T_up)
end
to_graphviz(boundary_conditions)

buoyancy_composition_diagram = @relation () begin
  momentum(V, v, b)

  # "Both T and S obey the tracer conservation equation"
  temperature(V, v, T, T_up)
  salinity(V, v, S, S_up)

  # "Buoyancy is determined from a linear equation of state"
  eos(b, T, S)

  bcs(S, T, S_up, T_up)
end
to_graphviz(buoyancy_composition_diagram, box_labels=:name, junction_labels=:variable, prog="fdp", graph_attrs=Dict(["sep" => "1.5"]))

buoyancy_cospan = oapply(buoyancy_composition_diagram, [
  Open(momentum,          [:V, :v, :b]),
  Open(tracer,            [:V, :v, :c, :c_up]),
  Open(tracer,            [:V, :v, :c, :c_up]),
  Open(equation_of_state, [:b, :T, :S]),
  Open(boundary_conditions, [:S, :T, :Ṡ, :T_up])])

buoyancy = apex(buoyancy_cospan)
to_graphviz(buoyancy)

buoyancy = expand_operators(buoyancy)
infer_types!(buoyancy)
resolve_overloads!(buoyancy)
to_graphviz(buoyancy)

s′ = loadmesh(Rectangle_30x10())
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())

include("../../grid_meshes.jl")
s′ = triangulated_grid(80,80, 10, 10, Point3D)
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())
xmax = maximum(x -> x[1], point(s))
zmax = maximum(x -> x[2], point(s))
wireframe(s)

# Some extra operators that haven't been up-streamed into CombinatorialSpaces yet.
include("dec_operators.jl")
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :L₁ => (x,y) -> L₁′(x,y,sd,hodge)
    :L₂ᵈ => (x,y) -> lie_derivative_flat(2, sd, x, y)
    :force => x -> x
    :∂_spatial => x -> begin
      left = findall(x -> x[1] ≈ 0.0, point(s))
      right = findall(x -> x[1] ≈ xmax, point(s))
      x[left] .= 276
      x[right] .= 279
      x
    end
    _ => default_dec_generate(sd, my_symbol, hodge)
  end
  return (args...) -> op(args...)
end

sim = eval(gensim(buoyancy))
fₘ = sim(s, generate)

S = map(point(s)) do (_,_,_)
  35.0
end
T = map(point(s)) do (x,z,_)
  273.15 + 4 + ((zmax-z)^2 + (xmax-x)^2)^(1/2)/(1e2)
end
extrema(T)
mesh(s′, color=T, colormap=:jet)
p = zeros(nv(s))
f = zeros(nv(s))
Fₛ = zeros(nv(s))
Fₜ = zeros(nv(s))
f = zeros(nv(s))
Cₛ = zeros(nv(s))
Cₜ = zeros(nv(s))

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
  VectorForm(Cₜ),
  VectorForm(S),
  VectorForm(Fₛ),
  VectorForm(qₛ),
  VectorForm(Cₛ)], Float64[], [
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
  :temperature_C,
  :S,
  :salinity_F,
  :salinity_q,
  :salinity_C])

gᶜ = 9.81
α = 2e-3
β = 5e-4
constants_and_parameters = (
  eos_g = gᶜ,
  eos_α = α,
  eos_β = β,
  momentum_U = t -> 0)

tₑ = 1e9

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
#prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
#prob = ODEProblem(fₘ, u₀, (0, 1e-3), constants_and_parameters)
#prob = ODEProblem(fₘ, u₀, (0, 1e-1), constants_and_parameters)
prob = ODEProblem(fₘ, u₀, (0, 1e2), constants_and_parameters)
soln = solve(prob, Tsit5(), dtmin=1e-3, force_dtmin=true)
soln.retcode != :Unstable || error("Solver was not stable")

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

mesh(s′, color=findnode(soln(2.2), :T), colormap=:jet)
extrema(findnode(soln(2.2), :T))

# Create a gif
begin
  frames = 100
  #fig, ax, ob = GLMakie.mesh(s′, color=findnode(soln(0), :T), colormap=:jet, colorrange=extrema(findnode(soln(tₑ), :h)))
  fig, ax, ob = GLMakie.mesh(s′, color=findnode(soln(0), :T), colormap=:jet, colorrange=extrema(findnode(soln(2.2), :T)))
  Colorbar(fig[1,2], ob)
  #record(fig, "oceananigans.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
  record(fig, "oceananigans.gif", range(0.0, 2.2; length=frames); framerate = 30) do t
    ob.color = findnode(soln(t), :T)
  end
end

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
