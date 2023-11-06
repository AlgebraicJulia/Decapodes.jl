using MLStyle
using Decapodes
using Decapodes.decapodes

function pprint(io::IO, exp::DecaExpr, pad=0) 
  pprint(io, "Context:\n", pad)
  pprint(io, exp.context, pad+2)
  pprint(io, "Equations:\n", pad)
  pprint(io, exp.equations, pad+2)
end

pprint(io::IO, exp::AbstractString, pad=0) = print(io, " "^pad*exp)

pprint(io::IO, exp::AbstractVector, pad=0) = map(exp) do x
  pprint(io, x, pad)
  println(io, "")
end

pprint(io::IO, exp::Equation, pad=0) = begin
  pprint(io, exp.lhs, pad)
  pprint(io, " = ", pad)
  pprint(io, exp.rhs, pad)
end

pprint(io::IO, exp::Judgement, pad=0) = begin
  pprint(io, "$(exp.var)::$(exp.dim) over $(exp.space)", pad)
end
pprint(io::IO, exp::Term, pad=0) = begin
  let ! = x -> sprint((io, y)->pprint(io, y, 0), x)
  @match exp begin
    Var(name) => print(io, name)
    Lit(name) => print(io, name)
    AppCirc1(fs, arg) => print(io, "($fs)($(!arg))")
    App1(f, arg) => print(io, "$f($(!arg))")
    App2(:*, arg1, arg2) => print(io, "$(!arg1) * $(!arg2)")
    App2(:+, arg1, arg2) => print(io, "$(!arg1) + $(!arg2)")
    App2(:-, arg1, arg2) => print(io, "$(!arg1) - $(!arg2)")
    App2(:/, arg1, arg2) => print(io, "$(!arg1) / $(!arg2)")
    App2(:./, arg1, arg2) => print(io, "$(!arg1) ./ $(!arg2)")
    App2(f, arg1, arg2) => print(io, "$f($(!arg1), $(!arg2))")
    Plus(args) => print(io, "$(join(map(!, args), " + "))")
    Mult(args) => print(io, "($(join(map(!, args), " * "))")
    Tan(var) => print(io, "∂ₜ($(!var))")
    _ => error("printing $exp")
  end
end
end