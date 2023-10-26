
abstract type AbstractDecapodeMorphism end

struct BCMorphism <: AbstractDecapodeMorphism
  morphism::ACSetTransformation
end

struct ICMorphism <: AbstractDecapodeMorphism
  morphism::ACSetTransformation
end

abstract type AbstractCollage end

struct Collage <: AbstractCollage
  bc::BCMorphism
  ic::ICMorphism
end

"""    function collate(dm::BCMorphism)

"Compile" a collage of Decapodes to a simulatable one.
```
"""
function collate(dm::BCMorphism)
  dm = dm.morphism
  d = SummationDecapode{Any, Any, Symbol}()
  copy_parts!(d, dm.codom, (:Var, :TVar, :Op1, :Op2, :Σ, :Summand))

  for (i,x) in enumerate(dm.components.Var.func)
    op_name = Symbol("∂_mask")
    mask_var = add_part!(d, :Var, type=dm.dom[i, :type], name=dm.dom[i, :name])
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
