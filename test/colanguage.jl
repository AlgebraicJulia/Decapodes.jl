# Test Decapode -> DecaExpr conversion on an empty Decapode.
empty_terms = parse_decapode(quote end)
EmptyDecapode = SummationDecapode(empty_terms)
@test empty_terms == Term(EmptyDecapode)

# Test Decapode -> DecaExpr conversion on a so-called discrete Decapode.
discrete_terms = parse_decapode(quote
  A::Form0
  B::Form1
  C::Form2
  (D,E)::Form0
end)
DiscreteDecapode = SummationDecapode(discrete_terms)
# Note that the roundtripping of the Judgements need not preserve order.
@test Term(DiscreteDecapode) == discrete_terms

# Test Decapode -> DecaExpr conversion on Op1s.
op1_terms = parse_decapode(quote
  (A,B)::Form0
  A == f(B)
end)
Op1Decapode = SummationDecapode(op1_terms)
@test Term(Op1Decapode) == op1_terms

# Test Decapode -> DecaExpr conversion on Op2s.
op2_terms = parse_decapode(quote
  (A,B,C)::Form0
  A == g(B,C)
end)
Op2Decapode = SummationDecapode(op2_terms)
@test Term(Op2Decapode) == op2_terms

# Test Decapode -> DecaExpr conversion on TVars.
dt_terms = parse_decapode(quote
  (A,B)::Form0
  A == ∂ₜ(B)
end)
dtDecapode = SummationDecapode(dt_terms)
@test Term(dtDecapode) == dt_terms

# Test Decapode <-> DecaExpr roundtripping on TVars.
@test dtDecapode == SummationDecapode(Term(dtDecapode))

# Test Decapode -> DecaExpr conversion on implicit TVars.
implicit_terms = parse_decapode(quote
  (A,B,C)::Form0
  ∂ₜ(A) == g(B,C)
end)
ImplicitDecapode = SummationDecapode(implicit_terms)
# This equality test fails, but note that they roundtrip to the same Decapode:
#@test Term(ImplicitDecapode) == implicit_terms
@test SummationDecapode(Term(ImplicitDecapode)) == ImplicitDecapode


# Test: Roundtripping James' case:
mixed_exp = parse_decapode(quote
  A::Form0{X}
  B::Form1{X}
  C::Form0{X}
  D::Form0{X}

  B == grad(A)
  C == f(A,B)
  ∂ₜ(A) == C
  ∂ₜ(D) == C + D
end)

Mix = SummationDecapode(mixed_exp)

@test SummationDecapode(Term(Mix)) == Mix

# Test Decapode -> DecaExpr conversion on Halfar.
halfar_terms = parse_decapode(quote
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯(d(h)))^(n-1)) * avg₀₁(h^(n+2)))
end)
HalfarDecapode = SummationDecapode(halfar_terms)
# Observe that expand_operators is called:
@test SummationDecapode(Term(HalfarDecapode)) == expand_operators(HalfarDecapode)
