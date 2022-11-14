using AlgebraicRewriting
using AlgebraicRewriting: Var as ARVar

@present SchDecapode(FreeSchema) begin
    (Var, TVar, Op1, Op2)::Ob
    (Type, Operator)::AttrType
    src::Hom(Op1, Var)
    tgt::Hom(Op1, Var)
    proj1::Hom(Op2, Var)
    proj2::Hom(Op2, Var)
    res::Hom(Op2, Var)
    incl::Hom(TVar, Var)
    
    op1::Attr(Op1, Operator)
    op2::Attr(Op2, Operator)
    type::Attr(Var, Type)
end

@present SchNamedDecapode <: SchDecapode begin
    Name::AttrType
    name::Attr(Var, Name)
end

@abstract_acset_type AbstractDecapode
@abstract_acset_type AbstractNamedDecapode <: AbstractDecapode

@acset_type Decapode(SchDecapode,
  index=[:src, :tgt, :res, :incl, :op1, :op2, :type]) <: AbstractDecapode

@acset_type NamedDecapode(SchNamedDecapode,
  index=[:src, :tgt, :res, :incl, :op1, :op2, :type, :name]) <: AbstractNamedDecapode

"""    fill_names!

add new variable names to all the variables that don't have names.
"""
function fill_names!(d::AbstractNamedDecapode)
    bulletcount = 1
    for i in parts(d, :Var)
        if !isassigned(d[:,:name],i)
            d[i,:name] = Symbol("•$bulletcount")
            bulletcount += 1
        end
    end
    for e in incident(d, :∂ₜ, :op1)
        s = d[e,:src]
        t = d[e, :tgt]
        d[t, :name] = append_dot(d[s,:name])
    end
    return d
end

function make_sum_unique!(d::AbstractNamedDecapode)
  num = 1
  for (i, name) in enumerate(d[:name])
    if(name == :sum)
      d[i, :name] = Symbol(join([String(name), string(num)] , "_"))
      num += 1
    end
  end
end

function expand_operators(d::AbstractNamedDecapode)
  e = SummationDecapode{Symbol, Symbol, Symbol}()
  copy_parts!(e, d, (:Var, :TVar, :Op2))
  expand_operators!(e, d)
  return e
end


function expand_operators!(e::AbstractNamedDecapode, d::AbstractNamedDecapode)
  newvar = 0
  for op in parts(d, :Op1)
    if !isa(d[op,:op1], AbstractArray)
      add_part!(e, :Op1, op1=d[op,:op1], src=d[op, :src], tgt=d[op,:tgt])
    else
      for (i, step) in enumerate(d[op, :op1])
        if i == 1
          newvar = add_part!(e, :Var, type=:infer, name=Symbol("•_$(op)_$(i)"))
          add_part!(e, :Op1, op1=step, src=d[op, :src], tgt=newvar)
        elseif i == length(d[op, :op1])
          add_part!(e, :Op1, op1=step, src=newvar, tgt=d[op,:tgt])
        else
          newvar′ = add_part!(e, :Var, type=:infer, name=Symbol("•_$(op)_$(i)"))
          add_part!(e, :Op1, op1=step, src=newvar, tgt=newvar′)
          newvar = newvar′
        end
      end
    end
  end
  return newvar
end
@present SchSummationDecapode <: SchNamedDecapode begin
  # Σ are the white nodes in the Decapode drawing
  # Summands are the edges that connect white nodes to variables (the projection maps)
  # because addition is commutative, we don't need to distinguish the order
  (Σ, Summand)::Ob
  summand::Hom(Summand, Var)
  summation::Hom(Summand, Σ)
  sum::Hom(Σ, Var)
end

@acset_type SummationDecapode(SchSummationDecapode,
  index=[:src, :tgt, :res, :incl, :op1, :op2, :type]) <: AbstractNamedDecapode


function expand_operators(d::SummationDecapode)
  e = SummationDecapode{Symbol, Symbol, Symbol}()
  copy_parts!(e, d, (:Var, :TVar, :Op2, :Σ, :Summand))
  expand_operators!(e, d)
  return e
end

function add_constant!(d::AbstractNamedDecapode, k::Symbol)
    return add_part!(d, :Var, type=:Constant, name=k)
end

function add_parameter(d::AbstractNamedDecapode, k::Symbol)
    return add_part!(d, :Var, type=:Parameter, name=k)
end


"""
These are the default rewrite rules used to do type inference.
"""
default_type_inference_rules = [
  begin
    L = @acset SummationDecapode{Any, Any, Any} begin
      Var = 2
      TVar = 1
      Op1 = 1
      
      type = [ARVar(:src_form), :infer]
      name = [ARVar(:src_name), ARVar(:tgt_name)]
      incl = [2]
      src = [1]
      tgt = [2]
      op1 = [:∂ₜ]
    end
    I = @acset SummationDecapode{Any, Any, Any} begin
      Var = 1
      
      type = [ARVar(:src_form)]
      name = [ARVar(:src_name)]
    end
    R = @acset SummationDecapode{Any, Any, Any} begin
      Var = 2
      TVar = 1
      Op1 = 1
      
      type = [ARVar(:src_form), ARVar(:src_form)]
      name = [ARVar(:src_name), ARVar(:tgt_name)]
      incl = [2]
      src = [1]
      tgt = [2]
      op1 = [:∂ₜ]
    end
    hil = AlgebraicRewriting.homomorphism(I, L)
    hir = AlgebraicRewriting.homomorphism(I, R)
    Rule(hil, hir)
  end]

function infer_types(d::SummationDecapode, rules::Vector{Rule{:DPO}})
  # Step 1: Convert to {Any,Any,Any} so we can use AlgebraicRewriting's Var.
  #d′ = migrate(SummationDecapode{Any, Any, Any}, d,
  #  Dict(:Var => :Var, :TVar => :TVar, :Op1 => :Op1, :Op2 => :Op2, :Σ => :Σ, :Summand => :Summand, :Type => :Type, :Operator => :Operator, :Name => :Name),
  #  Dict(:src => :src, :tgt => :tgt, :proj1 => :proj1, :proj2 => :proj2, :res => :res, :incl => :incl, :op1 => :op1, :op2 => :op2, :type => :type, :name => :name, :summand => :summand, :summation => :summation, :sum => :sum))
  d′ = migrate(SummationDecapode{Any, Any, Any}, d,
    merge(
      Dict(SchSummationDecapode.generators.Ob .=> SchSummationDecapode.generators.Ob),
      Dict(SchSummationDecapode.generators.AttrType .=> SchSummationDecapode.generators.AttrType)),
    merge(
      Dict(SchSummationDecapode.generators.Hom .=> SchSummationDecapode.generators.Hom),
      Dict(SchSummationDecapode.generators.Attr .=> SchSummationDecapode.generators.Attr)))

  # Step 2: Apply rules
  seq = Schedule[]
  append!(seq, RuleSchedule.(rules))
  #seq = Vector{Schedule}([RuleSchedule.(rules)])
  ar_step = ListSchedule(seq)
  end_condition(prev, curr) = :infer ∉ curr[:type]
  overall = WhileSchedule(ar_step, :main, end_condition)
  trajectory = apply_schedule(overall, G=d′)
  res = last(trajectory).G
  # Step 3: Convert back to {Any,Any,Symbol}.
  #migrate(SummationDecapode{Any, Any, Symbol}, res
  #  Dict(:Var => :Var, :TVar => :TVar, :Op1 => :Op1, :Op2 => :Op2, :Σ => :Σ, :Summand => :Summand, :Type => :Type, :Operator => :Operator, :Name => :Name),
  #  Dict(:src => :src, :tgt => :tgt, :proj1 => :proj1, :proj2 => :proj2, :res => :res, :incl => :incl, :op1 => :op1, :op2 => :op2, :type => :type, :name => :name, :summand => :summand, :summation => :summation, :sum => :sum))
  migrate(SummationDecapode{Any, Any, Symbol}, res,
    merge(
      Dict(SchSummationDecapode.generators.Ob .=> SchSummationDecapode.generators.Ob),
      Dict(SchSummationDecapode.generators.AttrType .=> SchSummationDecapode.generators.AttrType)),
    merge(
      Dict(SchSummationDecapode.generators.Hom .=> SchSummationDecapode.generators.Hom),
      Dict(SchSummationDecapode.generators.Attr .=> SchSummationDecapode.generators.Attr)))
end

infer_types(d::SummationDecapode) =
  infer_types(d, default_type_inference_rules)

