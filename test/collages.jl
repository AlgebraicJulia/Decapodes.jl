using Test
using Decapodes
using Catlab

# TODO: Special names for morphisms that match a method of grabbing boundary
# simplices.

# TODO: Initial conditions.
# - This will mean that we need to return both a masked Decapode, and a means of pluggin initial data in for physical quantities, replacing `constuct`.
# TODO: Temporal boundary conditions.
# TODO: General boundaries i.e. arbitrary functions of solutions.

# TODO: Add test with empty boundary Decapode.

# Test simple boundary masks.
DiffusionDynamics = @decapode begin
  K::Form0
  ∂ₜ(K) == ∘(d,⋆,d,⋆)(K)
end
DiffusionBoundaries = @decapode begin
  (Kb1, Kb2, Null)::Form0
end

# Test that simple boundary masks work on state variables.
StateMorphism = ACSetTransformation(
  DiffusionBoundaries, DiffusionDynamics,
  Var = [1,1,1])

DiffusionCollage = Decapodes.collate(StateMorphism)

@test DiffusionCollage == @acset SummationDecapode{Any, Any, Symbol} begin
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
  op2 = [:∂bKb1, :∂bKb2, :∂bNull]
  type  = [:Form0, :infer, :Parameter, :Form0, :Parameter, :Form0, :Parameter, :Form0]
  name  = [:K1, :K̇, :Kb1, :K2, :Kb2, :K3, :Null, :K]
end

# Test that simple boundary masks work on state and tangent variables.
StateTangentMorphism = ACSetTransformation(
  DiffusionBoundaries, DiffusionDynamics,
  Var = [1,1,2])

DiffusionCollage = Decapodes.collate(StateTangentMorphism)

@test DiffusionCollage == @acset SummationDecapode{Any, Any, Symbol} begin
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
  op2  = [:∂bKb1, :∂bKb2, :∂bNull]
  type  = [:Form0, :infer, :Parameter, :Form0, :Parameter, :Form0, :Parameter, :infer]
  name = [:K1, :K̇3, :Kb1, :K2, :Kb2, :K, :Null, :K̇]
end

# Test gensim on a collage.
#c = Collage(DiffusionDynamics, DiffusionBoundaries,
#  DiffusionMorphism, DiffusionSymbols)
#
#@test gensim(c) == gensim(DiffusionCollage)

# TODO: Test intermediate variable masks.
