# Import Dependencies 

## AlgebraicJulia Dependencies
using ACSets
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes

## External Dependencies
using ComponentArrays
using GLMakie
using GeometryBasics: Point3
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using Random
Point3D = Point3{Float64}

Random.seed!(0)

# Define Models

## Define vorticity streamflow formulation

Eq11InviscidPoisson = @decapode begin
  d𝐮::DualForm2
  𝐮::DualForm1
  ψ::Form0

  ψ == Δ₀⁻¹(⋆(d𝐮))
  𝐮 == ⋆(d(ψ))

  ∂ₜ(d𝐮) ==  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end

to_graphviz(Eq11InviscidPoisson)

## Apply boundary conditions with a collage

VorticityBoundaries = @decapode begin
  U::DualForm1
end
VorticityMorphism = @relation () begin
  bound_walls(C, Cb)
end
VorticitySymbols = Dict(
  :C => :𝐮,
  :Cb => :U)
VorticityBounded = collate(
  Eq11InviscidPoisson,
  VorticityBoundaries,
  VorticityMorphism,
  VorticitySymbols)

to_graphviz(VorticityBounded)

## Define phase segmentation process

CahnHilliard = @decapode begin
  C::DualForm0
  𝐯::DualForm1
  ∂ₜ(C) == 0.5 * ∘(⋆₁⁻¹,d₁,⋆₂)(
    d̃₀(C^3 - C - 0.5 * Δᵈ₀(C)) +
    C ∧ᵈᵈ₀₁ 𝐯)
end

to_graphviz(CahnHilliard)

## Compose bounded Navier-Stokes with phase field

NSPhaseFieldDiagram = @relation () begin
  navierstokes(𝐮)

  phasefield(𝐮)
end

draw_composition(NSPhaseFieldDiagram)

vort_ch = apex(oapply(NSPhaseFieldDiagram,
  [Open(VorticityBounded, [:𝐮]),
   Open(CahnHilliard, [:𝐯])]))

to_graphviz(vort_ch)

vort_ch = (resolve_overloads! ∘ infer_types! ∘ expand_operators)(vort_ch)

to_graphviz(vort_ch)

# Define the mesh

s = triangulated_grid(1.0, 1.0, 0.2, 0.2, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

wireframe(s, linewidth=2)
wireframe!(sd)

# Define constants, parameters, and initial conditions

## This is a dual 2-form, with values at the dual cells around primal vertices.
★ = dec_hodge_star(0,sd)
d𝐮₀ = ★ * ones(nv(sd))

## This is a dual 1-form, with values orthogonal to primal edges.
U₀ = zeros(ne(sd))

## This is a dual 0-form, with values at the dual vertices at the centers of triangles.
C₀ = sign(2,sd) .* rand(ntriangles(sd))

## Store these values to be passed to the solver.
u₀ = ComponentArray(
  navierstokes_d𝐮 = d𝐮₀,
  navierstokes_U = U₀,
  phasefield_C=C₀)

constants_and_parameters = (nothing=nothing,)

# Define how symbols map to Julia functions

boundary_triangles = boundary_inds(Val{2}, sd)

function simple_bounds(form,bvals)
  form[boundary_triangles] = bvals[boundary_triangles]
  form
end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :bound_walls => simple_bounds
    x => error("$x not matched")
  end
  return (args...) -> op(args...)
end

# Generate simulation 

## Write the simulation code to a file.
open("collage_mpf.jl", "w") do f
  write(f, string(gensim(vort_ch)))
end
sim = include("collage_mpf.jl") # in VSCode: include("../../collage_mpf.jl")

## Generate the simulation
fₘ = sim(sd, generate)

## Run simulation 

tₑ = 1e1

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-9), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

@save "collage_mpf.jld2" soln

# Visualize 
★ = dec_inv_hodge_star(0,sd)
f = Figure()
ax = GLMakie.Axis(f[1,1])
sctr = scatter!(ax, point(sd), color= ★ * soln(tₑ).navierstokes_d𝐮)
Colorbar(f[1,2], sctr)
f

