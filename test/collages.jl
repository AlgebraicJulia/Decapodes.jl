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

DiffusionMorphism = @relation () begin
  rb1_leftwall(C, Cb1)
  rb2_rightwall(C, Cb2)
  rb3(Ċ, Zero)
end

DiffusionCollage = collate(
  DiffusionDynamics,
  DiffusionBoundaries,
  DiffusionMorphism,
  Dict(
    :C => :K,
    :Ċ => :K̇,
    :Cb1 => :Kb1,
    :Cb2 => :Kb2,
    :Zero => :Null))
#to_graphviz(DiffusionCollage)

@test DiffusionCollage == @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 8
  TVar = 1
  Op1 = 2
  Op2 = 3
  src  = [1, 1]
  tgt  = [2, 2]
  proj1  = [3, 5, 7]
  proj2  = [4, 6, 8]
  res  = [1, 3, 2]
  incl  = [2]
  op1  = Any[:∂ₜ, [:d, :⋆, :d, :⋆]]
  op2  = [:rb1_leftwall, :rb2_rightwall, :rb3]
  type  = [:Form0, :infer, :Form0, :Form0, :Form0, :Form0, :infer, :infer]
  name  = [:r_K, :r_K̇, :r_K, :Kb1, :K, :Kb2, :K̇, :Null]
end

# Note: Since the order does not matter in which rb1 and rb2 are applied, it
# seems informal to state that one goes before the other.
# It might be better to provide a semantics for incident edges a la:
#Diffusion = @decapode begin
#  C::Form0
#  ∂ₜ(C) == rb3(∘(d,⋆,d,⋆)(rb1(C)))
#  ∂ₜ(C) == rb3(∘(d,⋆,d,⋆)(rb2(C)))
#end
# Such a technique would preserve the technical definition of "collage".
