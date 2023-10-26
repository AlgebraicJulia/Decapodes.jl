using Test
using Decapodes
using Catlab
using CombinatorialSpaces
using GeometryBasics: Point2
Point2D = Point2{Float64}

# TODO: Test intermediate variable masks.
# TODO: Initial conditions.
# - This will mean that we need to return both a masked Decapode, and a means of pluggin initial data in for physical quantities, replacing `constuct`.
# TODO: Temporal boundary conditions.
# TODO: General boundaries i.e. arbitrary functions of solutions. (More compelx relationships between what the boundary "is" and how to interpret state to provide values there.)

# TODO: Add test with empty boundary Decapode.

# Test simple boundary masks.
DiffusionDynamics = @decapode begin
  K::Form0
  ∂ₜ(K) == ∘(d,⋆,d,⋆)(K)
end
DiffusionBoundaries = @decapode begin
  (Kb1, Kb2)::Constant
  Null::Parameter
end

# Test that simple boundary masks work on state variables.
StateMorphism = BCMorphism(ACSetTransformation(
  DiffusionBoundaries, DiffusionDynamics,
  Var = [1,1,1]))

Diffusion = Decapodes.collate(StateMorphism)

@test Diffusion == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 8
  TVar = 1
  Op1 = 2
  Op2 = 3
  src = [8, 1]
  tgt = [2, 2]
  proj1  = [4, 6, 8]
  proj2  = [3, 5, 7]
  res  = [1, 4, 6]
  incl = [2]
  op1 = Any[:∂ₜ, [:d, :⋆, :d, :⋆]]
  op2 = [:∂_mask, :∂_mask, :∂_mask]
  type  = [:Form0, :infer, :Constant, :Form0, :Constant, :Form0, :Parameter, :Form0]
  name  = [:K1, :K̇, :Kb1, :K2, :Kb2, :K3, :Null, :K]
end

# Test that simple boundary masks work on state and tangent variables.
StateTangentMorphism = BCMorphism(ACSetTransformation(
  DiffusionBoundaries, DiffusionDynamics,
  Var = [1,1,2]))

Diffusion = Decapodes.collate(StateTangentMorphism)

@test Diffusion == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 8
  TVar = 1
  Op1 = 2
  Op2 = 3
  src  = [6, 1]
  tgt  = [8, 2]
  proj1  = [4, 6, 2]
  proj2  = [3, 5, 7]
  res = [1, 4, 8]
  incl = [8]
  op1  = Any[:∂ₜ, [:d, :⋆, :d, :⋆]]
  op2  = [:∂_mask, :∂_mask, :∂_mask]
  type  = [:Form0, :infer, :Constant, :Form0, :Constant, :Form0, :Parameter, :infer]
  name = [:K1, :K̇3, :Kb1, :K2, :Kb2, :K, :Null, :K̇]
end

# Test gensim on a collage.
@test gensim(StateTangentMorphism) == gensim(Diffusion)

# Test collate with empty boundaries.
EmptyBoundaries = @decapode begin
end

EmptyMorphism = BCMorphism(ACSetTransformation(
  EmptyBoundaries, DiffusionDynamics))

Diffusion = Decapodes.collate(EmptyMorphism)
@test Diffusion == DiffusionDynamics

# Test initial conditions via collage.
DiffusionICs = @decapode begin
  (Kic)::Form0
end
SingleICMorphism = ICMorphism(ACSetTransformation(
  DiffusionICs, DiffusionDynamics,
  Var = [1]))

# Test the collage workflow with IC and BC morphisms, step-by-step.
DiffusionCollage = Collage(
  BCMorphism(ACSetTransformation(
    DiffusionBoundaries, DiffusionDynamics,
    Var = [1,1,2])),
  ICMorphism(ACSetTransformation(
    DiffusionICs, DiffusionDynamics,
    Var = [1])))

# TODO: Make changes to the function that takes `generate` that automatically generates the mask application functions `∂Kb1`.
# TODO: Actually, it sounds like this can now be a single function called
# `∂_mask`.
"""    function make_bc_loader(dm::Collage)

Given a collage, return a function that accepts a mapping of Vars to masking functions.

This returned function will take in a mesh and non-boundary constants-and-parameters, and return a named tuple containing the parameters as evaluated on the mesh, and the regular parameters. (This final named tuple is suitable to be passed to an ODEProblem.)
"""
function make_bc_loader(dm::Collage)
  function loader(mask_funcs::Dict{Symbol, N}) where {N <: Function}
    function generator(sd, cs_ps::NamedTuple)
      vars = keys(mask_funcs)
      mask_pairs = map(values(mask_funcs)) do func
        func(sd)
      end
      merge(
        NamedTuple{Tuple(vars)}(mask_pairs),
        cs_ps)
    end
  end
  loader
end

function make_ic_loader(dm::Collage)
end

"""    function simulation_helper(dm::BCMorphism)

Given a collage, return functions to load boundary conditions, initial conditions, and a simulation.
```
"""
simulation_helper(dm::Collage) = (
  make_bc_loader(dm.bc), make_ic_loader(dm.ic), gensim(dm.bc))

# Boundary Condition:
#   Set all points with y-coordinate above cos(t) to sin(t).
# e.g. The polar ice cap is a source term, which expands and contracts
# seasonally, and so does its albedo.
constant_function(sd) = 
  t -> (map(x -> x[2] > cos(t), sd[:point]),
        fill(sin(t), count(x -> x[2] > 0.9, sd[:point])))
# Boundary Condition:
#   Set all points with y-coordinate below -0.9 to 0.0.
time_varying_function(sd) = (
  map(x -> x[2] > 0.9, sd[:point]),
  fill(0.0, count(x -> x[2] > 0.9, sd[:point])))

bc_loader = make_bc_loader(DiffusionCollage)

bc_generator = bc_loader(Dict(
    :Kb1 => constant_function,
    :Kb2 => time_varying_function))

s = loadmesh(Icosphere(5))
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Barycenter())

p = bc_generator(sd, (a = 1.0, b = 2.0))
@test typeof(p) <:
  NamedTuple{(:Kb2, :Kb1, :a, :b), Tuple{
    # Constant: a bit mask and a vector of floats
    Tuple{Vector{Bool}, Vector{Float64}},
    # Parameter: a function that returns ''
    N,
    # "Regular" constants and parameters: which happen to be floats here
    Float64, Float64}} where {N <: Function, M <: Function}

# Test that the p functions are passed and evaluate as expected.
@test typeof(p.Kb1(π/4)) == Tuple{Vector{Bool}, Vector{Float64}}
# Some basic trigonometry to check the correct number of vertices are selected.
@test count(p.Kb1(π/4)[1]) / nv(sd) -
  ((1/2)*(sqrt(1-1^2)*1 + asin(1)) - 
    (1/2)*(sqrt(1-cos(π/4)^2)*cos(π/4) + asin(cos(π/4)))) < 1e-3
@test all(p.Kb1(π/4)[2] .== sin(π/4))

# Test that p does not "blow up" memory by making a copy of the mesh.
# This 1% mark is a heuristic.
@test sizeof(p) / Base.summarysize(p)  < 0.01
@test sizeof(p) / Base.summarysize(sd) < 0.01
