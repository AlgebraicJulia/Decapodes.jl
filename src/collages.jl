
#struct Collage
#  src::SummationDecapode{Any,Any,Symbol}
#  tgt::SummationDecapode{Any,Any,Symbol}
#  uwd::Catlab.Programs.RelationalPrograms.UntypedUnnamedRelationDiagram{Symbol, Symbol}
#  symbols::Dict{Symbol, Symbol}
#end
#
#collate(c::Collage) = collate(c.src, c.tgt, c.uwd, c.symbols)

"""    function collate(dm::ACSetTransformation)

"Compile" a collage of Decapodes to a simulatable one.
```
"""
function collate(dm::ACSetTransformation)
  d = SummationDecapode{Any, Any, Symbol}()
  copy_parts!(d, dm.codom, (:Var, :TVar, :Op1, :Op2, :Σ, :Summand))
  # TODO: Type the mask value as a Parameter.

  for (i,x) in enumerate(dm.components.Var.func)
    mask_name = dm.dom[i, :name]
    op_name = Symbol("∂b" * string(mask_name))
    mask_var = add_part!(d, :Var, type=:Parameter, name=dm.dom[i, :name])
    tgt_name = dm.codom[x, :name]
    tgt_idx = only(incident(d, tgt_name, :name))
    d[tgt_idx, :name] = Symbol(string(d[tgt_idx, :name]) * string(i))
    res_var = add_part!(d, :Var, type=dm.codom[x, :type], name=tgt_name)
    if isempty(incident(d, tgt_idx, :incl)) # State variable
      add_part!(d, :Op2, proj1=res_var, proj2=mask_var, res=tgt_idx, op2=op_name)

      tangent_op1s = filter(x -> d[x, :op1]==:∂ₜ, incident(d, tgt_idx, :src))
      isempty(tangent_op1s) && continue
      d[only(tangent_op1s), :src] = res_var
      d[incident(d, tgt_idx, :incl), :incl] = res_var
    else # Tangent variable
      add_part!(d, :Op2, proj1=tgt_idx, proj2=mask_var, res=res_var, op2=op_name)

      tangent_op1s = filter(x -> d[x, :op1]==:∂ₜ, incident(d, tgt_idx, :tgt))
      isempty(tangent_op1s) && continue
      d[only(tangent_op1s), :tgt] = res_var
      d[incident(d, tgt_idx, :incl), :incl] = res_var
    end
  end
  d
end
