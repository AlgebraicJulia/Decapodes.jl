using Test
using Decapodes
using Catlab
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.CategoricalAlgebra


# TODO: Special names for morphisms that match a method of grabbing boundary
# simplices.

# TODO: Initial conditions.
# - This will mean that we need to return both a masked Decapode, and a means of pluggin initial data in for physical quantities, replacing `constuct`.
# TODO: Temporal boundary conditions.
# TODO: General boundaries i.e. arbitrary functions of solutions.

# TODO: Add test with empty boundary Decapode.

# Test simple boundary masks.
DiffusionDynamics = @decapode begin
  C::Form0
  ∂ₜ(C) == ∘(d,⋆,d,⋆)(C)
end
DiffusionBoundaries = @decapode begin
  (Cb1, Cb2, Zero)::Form0
end

DiffusionMorphism = @relation () begin
  rb1(C, Cb1)
  rb2(C, Cb2)
  rb3(Ċ, Zero)
end

DiffusionCollage = collate(DiffusionMorphism,
  DiffusionDynamics, DiffusionBoundaries, [
  (:C, :Cb1),
  (:C, :Cb2),
  (:Ċ, :Zero)])

# Note: Since the order does not matter in which rb1 and rb2 are applied, it
# seems informal to state that one goes before the other.
# It might be better to provide a semantics for incident edges a la:
#Diffusion = @decapode begin
#  C::Form0
#  #∂ₜ(C) == rb3(∘(d,⋆,d,⋆)(rb1(C)))
#  #∂ₜ(C) == rb3(∘(d,⋆,d,⋆)(rb2(C)))
#end

Diffusion = @decapode begin
  C::Form0
  ∂ₜ(C) == rb3(∘(d,⋆,d,⋆)(rb2(rb1(C))))
end

@test_fails is_isomorphic(DiffusionCollage, Diffusion)
