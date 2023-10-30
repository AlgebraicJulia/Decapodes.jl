using Test
using Decapodes
using Catlab
using CombinatorialSpaces
using GeometryBasics: Point2
Point2D = Point2{Float64}
using MultiScaleArrays

# TODO: General boundaries i.e. arbitrary functions of solutions. (More complex
# relationships between what the boundary "is" and how to interpret state to
# provide values there.)

# Test `transfer_sources!` transfers summands, op1s, proj1s, and proj2s.
Test1 = @decapode begin
  To::Form0
  A == B+From+D
  X == f(From)
  Y == g(From,G)
  Z == h(F,From)
end
Decapodes.transfer_sources!(Test1,
  only(incident(Test1, :From, :name)), only(incident(Test1, :To, :name)))

@test Test1 == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 10
  Op1 = 1
  Op2 = 2
  Σ = 1
  Summand = 3
  src = [1]
  tgt = [7]
  proj1 = [1, 10]
  proj2 = [6, 1]
  res = [9, 8]
  incl = Int64[]
  summand = [3, 1, 5]
  summation = [1, 1, 1]
  sum = [2]
  op1 = [:f]
  op2 = [:g, :h]
  type = [:Form0, :infer, :infer, :infer, :infer, :infer, :infer, :infer, :infer, :infer]
  name = [:To, :A, :B, :From, :D, :G, :X, :Z, :Y, :F]
end

# Test `transfer_targets!` transfers sums, op1s, and op2s.
Test2 = @decapode begin
  To::Form0
  From == B+C+D
  From == f(E)
  From == g(F,G)
end
Decapodes.transfer_targets!(Test2,
  only(incident(Test2, :From, :name)), only(incident(Test2, :To, :name)))

@test Test2 == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 8
  Op1 = 1
  Op2 = 1
  Σ = 1
  Summand = 3
  src = [8]
  tgt = [1]
  proj1 = [7]
  proj2 = [6]
  res = [1]
  incl = Int64[]
  summand = [3, 4, 5]
  summation = [1, 1, 1]
  sum = [1]
  op1 = [:f]
  op2 = [:g]
  type = [:Form0, :infer, :infer, :infer, :infer, :infer, :infer, :infer]
  name = [:To, :From, :B, :C, :D, :G, :F, :E]
end

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

# Test that simple boundary masks work on second order tangent variables.
SecondDerivatives = @decapode begin
  (X,Y)::Form0
  ∂ₜ(X) == Z
  Z == Δ(X)
  ∂ₜ(Y) == Δ(Z)
end
SingleBC = @decapode begin
  M::Form0
end
to_graphviz(SecondDerivatives)
TangentTangentMorphism = BCMorphism(ACSetTransformation(
  SingleBC, SecondDerivatives,
  Var = [4]))

SecondDerivativesBounded = Decapodes.collate(TangentTangentMorphism)

@test SecondDerivativesBounded == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 6
  TVar = 2
  Op1 = 4
  Op2 = 1
  src = [1, 1, 2, 6]
  tgt = [6, 4, 3, 3]
  proj1 = [4]
  proj2 = [5]
  res = [6]
  incl = [6, 3]
  summand = Int64[]
  summation = Int64[]
  sum = Int64[]
  op1 = [:∂ₜ, :Δ, :∂ₜ, :Δ]
  op2 = [:∂_mask]
  type = [:Form0, :Form0, :infer, :infer, :Form0, :infer]
  name = [:X, :Y, :Ẏ, :Z1, :M, :Z]
end

# Test gensim on a collage.
@test gensim(StateTangentMorphism) == gensim(Diffusion)

# Test that simple boundary masks work on intermediate variables.
DiffusionExpanded = expand_operators(DiffusionDynamics)
IntermediateMorphism = BCMorphism(ACSetTransformation(
  DiffusionBoundaries, DiffusionExpanded,
  Var = [3,3,3]))

Diffusion = Decapodes.collate(IntermediateMorphism)
@test Diffusion == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 11
  TVar = 1
  Op1 = 5
  Op2 = 3
  src = [1, 1, 3, 4, 5]
  tgt = [2, 11, 4, 5, 2]
  proj1 = [7, 9, 11]
  proj2 = [6, 8, 10]
  res = [3, 7, 9]
  incl = [2]
  op1 = [:∂ₜ, :d, :⋆, :d, :⋆]
  op2 = [:∂_mask, :∂_mask, :∂_mask]
  type = [:Form0, :infer, :infer, :infer, :infer, :Constant, :infer, :Constant, :infer, :Parameter, :infer]
  name = [:K, :K̇, Symbol("•_2_11"), Symbol("•_2_2"), 
Symbol("•_2_3"), :Kb1, Symbol("•_2_12"), :Kb2, Symbol("•_2_13"), :Null, Symbol("•_2_1")]
end

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

# We will use this mesh to test that Boundary and Initial Condition generators
# generate appropriate data structures.
# The fifth subdivision is chosen so that numeric results of applying the
# boundary condition computations can be compared with analytic methods up to
# some number of digits.
s = loadmesh(Icosphere(5))
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Barycenter())

# Boundary Condition:
#   Set all points with y-coordinate above cos(t) to sin(t).
# e.g. The polar ice cap is a source term, which expands and contracts
# seasonally, and so does its albedo.
time_varying_bc_function(sd) = 
  t -> (map(x -> x[2] > cos(t), sd[:point]),
        fill(sin(t), count(x -> x[2] > 0.9, sd[:point])))
# Boundary Condition:
#   Set all points with y-coordinate below -0.9 to 0.0.
constant_bc_function(sd) = (
  map(x -> x[2] > 0.9, sd[:point]),
  fill(0.0, count(x -> x[2] > 0.9, sd[:point])))

bc_loader = Decapodes.make_bc_loader(DiffusionCollage, 2)

bc_generator = bc_loader(Dict(
    :Kb1 => constant_bc_function,
    :Kb2 => constant_bc_function,
    :Null => time_varying_bc_function))

# Test that the BC generator returns a tuple of the correct type.
p = bc_generator(sd, (a = 1.0, b = 2.0))
@test typeof(p) <:
  NamedTuple{(:Kb2, :Kb1, :Null, :a, :b), Tuple{
    # Constant: a bit mask and a vector of floats
    Tuple{Vector{Bool}, Vector{Float64}},
    Tuple{Vector{Bool}, Vector{Float64}},
    # Parameter: a function that returns ''
    N,
    # "Regular" constants and parameters: which happen to be floats here
    Float64, Float64}} where {N <: Function}

# Test that the p functions are passed and evaluate as expected.
@test typeof(p.Null(π/4)) == Tuple{Vector{Bool}, Vector{Float64}}
# Some basic trigonometry to check the correct number of vertices are selected.
@test count(p.Null(π/4)[1]) / nv(sd) -
  abs((1/2)*(sqrt(1-       1^2)*1        + asin(1       )) - 
      (1/2)*(sqrt(1-cos(π/4)^2)*cos(π/4) + asin(cos(π/4)))) < 1e-3
@test all(p.Null(π/4)[2] .== sin(π/4))

# Test that p does not "blow up" memory by making a copy of the mesh.
# This 1% mark is a heuristic.
@test sizeof(p) / Base.summarysize(p)  < 0.01
@test sizeof(p) / Base.summarysize(sd) < 0.01

# Test that the boundary condition generator will catch if initial conditions
# are not provided.
bc_generator = bc_loader(Dict(
    :Kb1 => constant_bc_function,
    #:Kb2 => constant_bc_function,
    :Null => time_varying_bc_function))
# Test that the BC generator returns a tuple of the correct type.
@test_throws "BC Variable Kb2 is not given a generating function." bc_generator(sd, (a = 1.0, b = 2.0))

# Initial Condition:
#   Assign each vertex the sin of their y-coordinate.
ic_function(sd) = map(x -> sin(x[2]), sd[:point])
# This function will return a vector of length the number of edges in the mesh,
# and we will erroneously assign it to a Form0 so as to test error-checking.
false_ic_function(sd) = ones(ne(sd))

ic_loader = Decapodes.make_ic_loader(DiffusionCollage, 2)
ic_generator = ic_loader(Dict(:Kic => ic_function))
u = ic_generator(sd)

# Test that the initial condition functions are passed and evaluate as expected.
@test typeof(u) == PhysicsState{VectorForm{Float64}, Float64}
@test findnode(u, :K) == ic_function(sd)
# Test that the generator is equivalent to setting initial conditions manually.
# Note: Testing equality with `==` fails for MultiScaleArrays.
# i.e. This test will fail:
#@test construct(PhysicsState,[VectorForm(ic_function(sd))], Float64[], [:K]) ==
#      construct(PhysicsState,[VectorForm(ic_function(sd))], Float64[], [:K])
# So, we check for equality between the accessed components.
v = construct(PhysicsState, [VectorForm(ic_function(sd))], Float64[], [:K])
@test all(findnode(u, :K) .- findnode(v, :K) .== 0)

# Test the error-checking on the length of initial condition data.
ic_generator = ic_loader(Dict(:Kic => false_ic_function))
@test_throws "IC Variable Kic was declared to be a Form0, but is of length 7680, not 2562." ic_generator(sd)

# Test the `simulation_helper` function, which wraps the boundary conditions
# loader, initial conditions loader, and the simulation generator.
sim = gensim(Decapodes.collate(DiffusionCollage.bc))

auto_bc_loader, auto_ic_loader, auto_sim =
  simulation_helper(DiffusionCollage, dimension=2)
@test auto_bc_loader == bc_loader
@test auto_ic_loader == ic_loader
@test auto_sim == sim
