using Decapodes
using Catlab
using Catlab.CategoricalAlgebra

function Term(s::SummationDecapode)
  s = expand_operators(s)
  judgements = map(parts(s,:Var)) do v
    var = s[v, :name]
    typ = s[v, :type]
    Judgement(var, typ, :I)
  end

  op1s = map(parts(s, :Op1)) do op
    x = Var(s[op, [:src, :name]])
    y = Var(s[op, [:tgt, :name]])
    f = s[op, :op1]
    if f == :∂ₜ
      Eq(y, Tan(x))
    else
      Eq(y, App1(f, x))
    end
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
