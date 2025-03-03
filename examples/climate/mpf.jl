# Import Dependencies 

## AlgebraicJulia Dependencies
using ACSets
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes

## External Dependencies
using ComponentArrays
using CairoMakie
using Distributions
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
  μ::Form0

  ψ == Δ₀⁻¹(⋆(d𝐮))
  𝐮 == ⋆(d(ψ))

  ∂ₜ(d𝐮) == μ * ∘(⋆, d, ⋆, d)(d𝐮) - ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
end

to_graphviz(Eq11InviscidPoisson)
(to_graphviz ∘ resolve_overloads! ∘ infer_types! ∘ expand_operators)(Eq11InviscidPoisson)

## Apply boundary conditions with a collage

VorticityBoundaries = @decapode begin
  U::DualForm1
  DU::DualForm2
end
VorticityMorphism = @relation () begin
  bound_dual1form(Flow, FlowBoundaryValues)
  bound_dual2form(Vorticity, VorticityBoundaryValues)
end
VorticitySymbols = Dict(
  :Flow => :𝐮,
  :FlowBoundaryValues => :U,
  :Vorticity => :d𝐮,
  :VorticityBoundaryValues => :DU)
VorticityBounded = collate(
  Eq11InviscidPoisson,
  VorticityBoundaries,
  VorticityMorphism,
  VorticitySymbols)

to_graphviz(VorticityBounded)

## Define phase segmentation process

CahnHilliard = @decapode begin
  C::Form0
  𝐯::DualForm1
  (D,γ)::Constant
  η::Constant
  F::Constant
  ∂ₜ(C) == D * ∘(⋆,d,⋆)(
    F * d(C^3 - C - γ * Δ(C)) +
    η * (C ∧ ♭♯(𝐯)))
end

to_graphviz(CahnHilliard)

## Define viscosity

SigmoidalViscosity = @decapode begin
  (C,μ)::Form0
  (L,k,J)::Constant
  μ == L / (1 + exp(-(k)*C)) + J
end

to_graphviz(SigmoidalViscosity)

## Compose bounded Navier-Stokes with phase field

NSPhaseFieldDiagram = @relation () begin
  navierstokes(𝐮,μ)

  viscosity(μ,C)

  phasefield(𝐮,C)
end

draw_composition(NSPhaseFieldDiagram)

vort_ch = apex(oapply(NSPhaseFieldDiagram,
  [Open(VorticityBounded, [:𝐮,:μ]),
   Open(SigmoidalViscosity, [:μ,:C]),
   Open(CahnHilliard, [:𝐯,:C])]))

to_graphviz(vort_ch)

### Infer types and resolve operators here for illustrative purposes.
### This code is called internally in `gensim`.
vort_ch = (resolve_overloads! ∘ infer_types! ∘ expand_operators)(vort_ch)

to_graphviz(vort_ch)

# Define the mesh

s = triangulated_grid(1.0, 1.0, 0.0125, 0.0125, Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter())

function plot_primal_dual_wireframes()
  f = Figure()
  ax = CairoMakie.Axis(f[1,1])
  wireframe!(ax, s, linewidth=2)
  wireframe!(ax, sd)
  f
end
plot_primal_dual_wireframes()

# Define constants, parameters, and initial conditions

## d𝐮₀ is a dual 2-form, with values at the dual cells around primal vertices.
★0 = dec_hodge_star(0,sd)
distribution = MvNormal([0.5, 0.5, 0.0], Diagonal([1/8, 1/8, 1e-9]))
d𝐮₀ = normalize(★0 * map(x -> pdf(distribution, x), point(sd)), 1)

DU₀ = zeros(nv(sd))

## U₀ is a dual 1-form, with values orthogonal to primal edges.
U₀ = zeros(ne(sd))

## C₀ is a primal 0-form, with values at primal vertices.
C₀ = map(point(sd)) do (x,y,z)
  x < 0.5 + sin((y)*10)/4 ? -1 : +1
end

## Store these values to be passed to the solver.
u₀ = ComponentArray(
  navierstokes_d𝐮 = d𝐮₀,
  navierstokes_U = U₀,
  navierstokes_DU = DU₀,
  C = C₀)

constants_and_parameters = (
  viscosity_L = 9e-3,
  viscosity_k = 6,
  viscosity_J = 1e-3,
  phasefield_F = 1e-1,
  phasefield_D = 5e-2,
  phasefield_γ = (1e-2)^2,
  phasefield_η = 1e12)

# Define how symbols map to Julia functions

boundary_edges = boundary_inds(Val{1}, sd)
boundary_vertices = boundary_inds(Val{0}, sd)

function simple_dual1form_bounds(form, bvals)
  form[boundary_edges] = bvals[boundary_edges]
  form
end
function simple_dual2form_bounds(form, bvals)
  form[boundary_vertices] = bvals[boundary_vertices]
  form
end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :bound_dual1form => simple_dual1form_bounds
    :bound_dual2form => simple_dual2form_bounds
    :exp => x -> exp.(x)
    x => error("$x not matched")
  end
  return (args...) -> op(args...)
end

# Generate simulation 

## Write the simulation code to a file.
open("collage_mpf.jl", "w") do f
  write(f, string(gensim(vort_ch)))
end
sim = include("../../collage_mpf.jl") # At the terminal, use: sim = include("collage_mpf.jl")

## Generate the simulation
fₘ = sim(sd, generate)

## Run simulation 

tₑ = 2e1

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
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
function plot_star_vorticities()
  ★ = dec_inv_hodge_star(0,sd)
  f = Figure(size = (400, 1200))
  ax = CairoMakie.Axis(f[1,1], aspect=1, title="★(d𝐮(0))")
  msh = mesh!(ax, s, color= ★ * soln(0).navierstokes_d𝐮)
  Colorbar(f[1,2], msh)
  ax2 = CairoMakie.Axis(f[2,1], aspect=1, title="★(d𝐮(tₑ))")
  msh2 = mesh!(ax2, s, color= ★ * soln(tₑ).navierstokes_d𝐮)
  Colorbar(f[2,2], msh2)
  ax3 = CairoMakie.Axis(f[3,1], aspect=1, title="★(d𝐮(tₑ) - d𝐮(0))")
  msh3 = mesh!(ax3, s, color= ★ * (soln(tₑ).navierstokes_d𝐮 - soln(0).navierstokes_d𝐮))
  Colorbar(f[3,2], msh3)
  f
end
f = plot_star_vorticities()
save("star_vorticities.png", f)

function plot_final_phasefield()
  f = Figure()
  ax = CairoMakie.Axis(f[1,1], title="C(0)")
  msh = mesh!(ax, s, color=soln(tₑ).C)
  Colorbar(f[1,2], msh)
  f
end
f = plot_final_phasefield()
save("final_phasefield.png", f)

function animate_star_vorticity(file_name, soln)
  ★ = dec_inv_hodge_star(0,sd)
  time = Observable(0.0)
  fig = Figure()
  Label(fig[1, 1, Top()], @lift("★(d𝐮) at $($time)"), padding = (0, 0, 5, 0))
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s,
    color=@lift(★ * soln($time).navierstokes_d𝐮))
  Colorbar(fig[1,2], msh)
  record(fig, file_name, range(0, tₑ; length=40); framerate = 10) do t
    time[] = t
  end
end
animate_star_vorticity("star_vorticity.mp4", soln)

function animate_phasefield(file_name, soln)
  time = Observable(0.0)
  fig = Figure()
  Label(fig[1, 1, Top()], @lift("C at $($time)"), padding = (0, 0, 5, 0))
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s,
    color=@lift(soln($time).C))
  Colorbar(fig[1,2], msh)
  record(fig, file_name, range(0, tₑ; length=40); framerate = 10) do t
    time[] = t
  end
end
animate_phasefield("phasefield.mp4", soln)
