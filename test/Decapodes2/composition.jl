using Test
using Decapodes
using Catlab
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.CategoricalAlgebra
using Catlab.CSetDataStructures

import Decapodes: OpenSummationDecapode, OpenPode, oapply, oapply_rename
# @testset "Composition" begin
# Simplest possible decapode relation.
TrivialExprBody = quote
  H::Form0{X}
end

trivalExpr = parse_decapode(TrivialExprBody)
Trivial = SummationDecapode(trivalExpr)
trivial_relation = @relation () begin
  trivial(H)
end

otrivial = OpenPode(Trivial, [:H])

trivial_comp = oapply(trivial_relation, [OpenPode(Trivial, [:H])])
apex(trivial_comp)

@test apex(oapply(trivial_relation, OpenPode(Trivial, [:H]))) == Trivial
@test apex(trivial_comp) == Trivial

# Multiple variables and an equation.
AdvectionExprBody =  quote
  C::Form0{X}
  V::Form1{X}
  ϕ::Form1{X}
  ϕ == ∧₀₁(C,V)
end
advExpr = parse_decapode(AdvectionExprBody)
Advection = SummationDecapode(advExpr)
adv_relation = @relation () begin
  advection(C,V,ϕ)
end
adv_comp = oapply(adv_relation, [OpenPode(Advection, [:C,:V,:ϕ])])
adv_comp_expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 3
  type = [:Form0, :Form1, :Form1]
  name = [:C, :V, :ϕ]
  Op2 = 1
  proj1 = [1]
  proj2 = [2]
  res = [3]
  op2 = [:∧₀₁]
end
@test apex(adv_comp) == adv_comp_expected

# This is the example from the "Overview" page in the docs.
DiffusionExprBody =  quote
  C::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  ∘(k, d₀)(C)
end
AdvectionExprBody = quote
  C::Form0{X}
  V::Form1{X}
  ϕ::Form1{X}

  ϕ == ∧₀₁(C,V)
end
SuperpositionExprBody = quote
  C::Form0{X}
  Ċ::Form0{X}
  ϕ::Form1{X}
  ϕ₁::Form1{X}
  ϕ₂::Form1{X}

  ϕ == ϕ₁ + ϕ₂
  Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
  ∂ₜ(C) == Ċ
end

difExpr = parse_decapode(DiffusionExprBody)
Diffusion = SummationDecapode(difExpr)

advExpr = parse_decapode(AdvectionExprBody)
Advection = SummationDecapode(advExpr)

supExpr = parse_decapode(SuperpositionExprBody)
#Superposition = SummationDecapode(supExpr)
Superposition = SummationDecapode(supExpr)

compose_diff_adv = @relation (C,V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

decapodes_vars = [
  OpenPode(Diffusion, [:C, :ϕ]),
  OpenPode(Advection, [:C, :ϕ, :V]),
  OpenPode(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])]
debugg = oapply_rename(compose_diff_adv, decapodes_vars)

dif_adv_sup = oapply(compose_diff_adv, decapodes_vars)

dif_adv_sup_expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 6
  type = [:Form0, :Form1, :Form1, :Form1, :infer, :Form1]
  name = [:C, :ϕ₁, :V, :ϕ₂, :superposition_Ċ, :ϕ]

  TVar = 1
  incl = [5]

  Op1 = 3
  src = [1,6,1]
  tgt = [2,5,5]
  op1 = [[:k, :d₀], [:⋆₀⁻¹, :dual_d₁, :⋆₁], :∂ₜ]

  Op2 = 1
  proj1 = [1]
  proj2 = [3]
  res = [4]
  op2 = [:∧₀₁]

  Σ = 1
  sum = [6]

  Summand = 2
  summand = [2,4]
  summation = [1,1]
end
@test apex(dif_adv_sup) == dif_adv_sup_expected

# Test some other permutation of the symbols yields the same decapode.
compose_diff_adv = @relation (C,V) begin
  diffusion(C, ϕ₁)
  #advection(C, ϕ₂, V)
  #superposition(ϕ₁, ϕ₂, ϕ, C)
  advection(V, ϕ₂, C)
  superposition(ϕ₁, ϕ₂, C, ϕ)
end

decapodes_vars = [
  OpenPode(Diffusion, [:C, :ϕ]),
  #OpenPode(Advection, [:C, :ϕ, :V]),
  #OpenPode(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])]
  OpenPode(Advection, [:V, :ϕ, :C]),
  OpenPode(Superposition, [:ϕ₁, :ϕ₂, :C, :ϕ])]
dif_adv_sup = oapply(compose_diff_adv, decapodes_vars)

@test apex(dif_adv_sup) == dif_adv_sup_expected

# Test that Op2s are properly de-duplicated.
AdvectionExprBody = quote
  C::Form0{X}
  Ċ::Form0{X}
  V::Form1{X}
  ϕ::Form1{X}
  ϕ == ∧₀₁(C,V)
  ∂ₜ(C) == Ċ
  ∂ₜ(C) == ⋆(ϕ)
end

advExpr = parse_decapode(AdvectionExprBody)

Advection = expand_operators(SummationDecapode(advExpr))

self_adv = @relation () begin
  advection₁(C,V,ϕ)
  advection₂(C,V,ϕ)
end

#adv_adv = [
# OpenPode(Advection, [:C,:V,:ϕ]),
# OpenPode(Advection, [:C,:V,:ϕ])]
#adv_adv_comp = oapply(self_adv, adv_adv)
## De-duplicate Op1s.
#unique_by!(adv_adv_comp, :Op1, [:src, :tgt, :op1])
## De-duplicate Op2s.
#unique_by!(adv_adv_comp, :Op2, [:proj1, :proj2, :res, :op2])
#adv_adv_comp_expected = @acset SummationDecapode{Any, Any, Symbol} begin
#  Var = 3
#  type = [:Form0, :Form1, :Form1]
#  name = [:C, :V, :ϕ]
#  Op2 = 1
#  proj1 = [1]
#  proj2 = [2]
#  res = [3]
#  op2 = [:∧₀₁]
#end
#@test apex(adv_adv_comp) == adv_adv_comp_expected

# end
