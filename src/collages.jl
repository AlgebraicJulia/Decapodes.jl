
"""    function collate(acset, column_names::Vector{Symbol})

Create a collage of two Decapodes that simulates with boundary conditions.
```
"""
function collate(equations, boundaries, uwd, symbols)
  # TODO: This is assuming only "restriction"-type morphisms.

  f = SummationDecapode{Any, Any, Symbol}()
  # TODO: Double-check
  copy_parts!(f, equations, (:Var, :TVar, :Op1, :Op2, :Σ, :Summand))

  # TODO: Throw an error if the user tries to use a boundary value differential
  # form that is of a different type of the thing that we are applying the bound
  # to. i.e. Form1 but target is a Form0.

  # TODO: This sets restrictions as Op1s. They are actually Op2s. i.e. Use `bv`.
  for b in boxes(uwd)
    ps = incident(uwd, b, :box)
    ev = first(ps)
    bv = last(ps)
    en_key = uwd[junction(uwd, ev), :variable]
    bn_key = uwd[junction(uwd, bv), :variable]
    en = symbols[en_key]
    bn = symbols[bn_key]
    var = only(incident(f, en, :name))
    b_var = add_part!(f, :Var, type=f[var, :type], name=f[var, :name])
    f[var, :name] = Symbol("r$(b)_" * string(f[var, :name]))
    s_var = add_part!(f, :Var, type=boundaries[only(incident(boundaries, bn, :name)), :type], name=bn)
    add_part!(f, :Op2, proj1=b_var, proj2=s_var, res=var, op2=uwd[b, :name])

    # Update tangent variable pointers, if any.
    tangent_op1s = filter(x -> f[x, :op1]==:∂ₜ, incident(f, var, :src))
    isempty(tangent_op1s) && continue
    f[only(tangent_op1s), :src] = b_var
  end

  f
end
