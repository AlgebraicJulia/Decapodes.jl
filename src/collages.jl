
"""    function collate(acset, column_names::Vector{Symbol})

Create a collage of two Decapodes that simulates with boundary conditions.
```
"""
function collate(equations, boundaries, uwd, symbols)
  # TODO: This is assuming only "restriction"-type morphisms.

  f = SummationDecapode{Any, Any, Symbol}()
  # TODO: Double-check
  copy_parts!(f, equations, (:Var, :TVar, :Op1, :Op2, :Σ, :Summand, :Sum))

  # TODO: This sets restrictions as Op1s. They are actually Op2s. i.e. Use `bv`.
  for b in boxes(uwd)
    ps = incident(uwd, b :box)
    ev = first(ps)
    bv = last(ps)
    var = incident(f, ev, :name)
    b_var = add_part!(f, :Var, type=f[var, :type], name=f[var, :name])
    f[var, :name] = Symbol("r_" * string(f[var, :name]))
    add_part!(f, :Op1, src=b_var, tgt=var, op1=uwd[b, :name])
  end

  f
end
