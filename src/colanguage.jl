using Decapodes
using Catlab
using Catlab.CategoricalAlgebra

function Term(s::SummationDecapode)
  judgements = map(parts(s,:Var)) do v
    var = s[v, :name]
    typ = s[v, :type]
    Judgement(Var(var), typ, :X)
  end

  op1s = map(parts(s, :Op1)) do op
    x = Var(s[op, [:src, :name]])
    y = Var(s[op, [:tgt, :name]])
    f = s[op, :op1]
    if f == :∂ₜ
      y = Tan(y)
    end
    Eq(y, App1(f, x))
  end

  op2s = map(parts(s, :Op2)) do op
    x = Var(s[op, [:proj1, :name]])
    y = Var(s[op, [:proj2, :name]])
    z = Var(s[op, [:res, :name]])
    f = s[op, :op2]
    Eq(z, App2(f, x, y))
  end

  sums = map(parts(s, :Σ)) do σ
    terms = map(Var, s[incident(s, σ, :summation), [:summand, :name]])
    Eq(Var(s[σ, [:sum,:name]]), Plus(terms))
  end
  Decapodes.DecaExpr(judgements, vcat(op1s, op2s, sums))
end

dexp = parse_decapode(quote
  A::Form0{X}
  B::Form1{X}
  C::Form0{X}
  D::Form0{X}

  B == grad(A)
  C == f(A,B)
  ∂ₜ(A) == C
  ∂ₜ(D) == C + D
end)

d = SummationDecapode(dexp)

dexpr′ = Term(d)
d′ = SummationDecapode(Term(d))