
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
    bn = symbols[bn_key] # This will be used when we do error-checking on types.
    println(en_key, en)
    var = only(incident(f, en, :name))
    b_var = add_part!(f, :Var, type=f[var, :type], name=f[var, :name])
    f[var, :name] = Symbol("r_" * string(f[var, :name]))
    #add_part!(f, :Op1, src=b_var, tgt=var, op1=uwd[b, :name])
    s_var = add_part!(f, :Var, type=f[var, :type], name=bn)
    add_part!(f, :Op2, proj1=b_var, proj2=s_var, res=var, op2=uwd[b, :name])
  end

  f
end
