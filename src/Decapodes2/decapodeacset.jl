using DataStructures

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
    if !isassigned(d[:,:name],i) || isnothing(d[i, :name])
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

function make_sum_mult_unique!(d::AbstractNamedDecapode)
  snum = 1
  mnum = 1
  for (i, name) in enumerate(d[:name])
    if(name == :sum)
      d[i, :name] = Symbol("sum_$(snum)")
      snum += 1
    elseif(name == :mult)
      d[i, :name] = Symbol("mult_$(mnum)")
      mnum += 1
    end
  end
end

# Note: This hard-bakes in Form0 through Form2, and higher Forms are not
# allowed.
function recognize_types(d::AbstractNamedDecapode)
  unrecognized_types = setdiff(d[:type], [:Form0, :Form1, :Form2, :DualForm0,
                          :DualForm1, :DualForm2, :Literal, :Parameter,
                          :Constant, :infer])
  isempty(unrecognized_types) ||
    error("Types $unrecognized_types are not recognized.")
end

function expand_operators(d::AbstractNamedDecapode)
  #e = SummationDecapode{Symbol, Symbol, Symbol}()
  e = SummationDecapode{Any, Any, Symbol}()
  copy_parts!(e, d, (:Var, :TVar, :Op2))
  expand_operators!(e, d)
  return e
end


function expand_operators!(e::AbstractNamedDecapode, d::AbstractNamedDecapode)
  newvar = 0
  for op in parts(d, :Op1)
    if !isa(d[op,:op1], AbstractArray)
      add_part!(e, :Op1, op1=d[op,:op1], src=d[op, :src], tgt=d[op,:tgt])
    elseif length(d[op, :op1]) == 1
      add_part!(e, :Op1, op1=only(d[op,:op1]), src=d[op, :src], tgt=d[op,:tgt])
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
  #e = SummationDecapode{Symbol, Symbol, Symbol}()
  e = SummationDecapode{Any, Any, Symbol}()
  copy_parts!(e, d, (:Var, :TVar, :Op2, :Σ, :Summand))
  expand_operators!(e, d)
  return e
end

"""
  function contract_operators(d::SummationDecapode)

Find chains of Op1s in the given Decapode, and replace them with
a single Op1 with a vector of function names. After this process,
all Vars that are not a part of any computation are removed.
"""
function contract_operators(d::SummationDecapode)
  e = expand_operators(d)
  contract_operators!(e)
  #return e
end

function contract_operators!(d::SummationDecapode)
  chains = find_chains(d)
  filter!(x -> length(x) != 1, chains)
  for chain in chains
    add_part!(d, :Op1, src=d[:src][first(chain)], tgt=d[:tgt][last(chain)], op1=Vector{Symbol}(d[:op1][chain]))
  end
  rem_parts!(d, :Op1, sort!(vcat(chains...)))
  remove_neighborless_vars!(d)
end

"""
  function remove_neighborless_vars!(d::SummationDecapode)

Remove all Vars from the given Decapode that are not part of any computation.
"""
function remove_neighborless_vars!(d::SummationDecapode)
  neighborless_vars = setdiff(parts(d,:Var),
                              union(d[:src],
                                    d[:tgt],
                                    d[:proj1],
                                    d[:proj2],
                                    d[:res],
                                    d[:sum],
                                    d[:summand],
                                    d[:incl]))
  #union(map(x -> t5_orig[x], [:src, :tgt])...) alternate syntax
  #rem_parts!(d, :Var, neighborless_vars)
  rem_parts!(d, :Var, sort!(neighborless_vars))
  d
end

"""
  function find_chains(d::SummationDecapode)

Find chains of Op1s in the given Decapode. A chain ends when the
target of the last Op1 is part of an Op2 or sum, or is a target
of multiple Op1s.
"""
function find_chains(d::SummationDecapode)
  chains = []
  visited = falses(nparts(d, :Op1))
  # TODO: Re-write this without two reduce-vcats.
  chain_starts = unique(reduce(vcat, reduce(vcat,
                        #[incident(d, Decapodes.infer_states(d), :src),
                        [incident(d, Vector{Int64}(filter(i -> !isnothing(i), Decapodes.infer_states(d))), :src),
                         incident(d, d[:res], :src),
                         incident(d, d[:sum], :src)])))
  
  s = Stack{Int64}()
  foreach(x -> push!(s, x), chain_starts)
  while !isempty(s)
    # Start a new chain.
    op_to_visit = pop!(s)
    curr_chain = []
    while true
      visited[op_to_visit] = true
      append!(curr_chain, op_to_visit)

      tgt = d[op_to_visit, :tgt]
      next_op1s = incident(d, tgt, :src)
      next_op2s = vcat(incident(d, tgt, :proj1), incident(d, tgt, :proj2))
      if (length(next_op1s) != 1 ||
          length(next_op2s) != 0 ||
          is_tgt_of_many_ops(d, tgt) ||
          !isempty(incident(d, tgt, :sum)) ||
          !isempty(incident(d, tgt, :summand)))
        # Terminate chain.
        append!(chains, [curr_chain])
        for next_op1 in next_op1s
          visited[next_op1] || push!(s, next_op1)
        end
        break
      end
      # Try to continue chain.
      op_to_visit = only(next_op1s)
    end
  end
  return chains
end

function add_constant!(d::AbstractNamedDecapode, k::Symbol)
    return add_part!(d, :Var, type=:Constant, name=k)
end

function add_parameter(d::AbstractNamedDecapode, k::Symbol)
    return add_part!(d, :Var, type=:Parameter, name=k)
end

# TODO: You could write a method which auto-generates these rules given degree N.
"""
These are the default rules used to do type inference in the 1D exterior calculus.
"""
op1_inf_rules_1D = [
  # Rules for ∂ₜ 
  (src_type = :Form0, tgt_type = :Form0, op_names = [:∂ₜ]),
  (src_type = :Form1, tgt_type = :Form1, op_names = [:∂ₜ]),

  # Rules for d
  (src_type = :Form0, tgt_type = :Form1, op_names = [:d, :d₀]),
  (src_type = :DualForm0, tgt_type = :DualForm1, op_names = [:d, :dual_d₀, :d̃₀]),

  # Rules for ⋆
  (src_type = :Form0, tgt_type = :DualForm1, op_names = [:⋆, :⋆₀]),
  (src_type = :Form1, tgt_type = :DualForm0, op_names = [:⋆, :⋆₁]),
  (src_type = :DualForm1, tgt_type = :Form0, op_names = [:⋆, :⋆₀⁻¹]),
  (src_type = :DualForm0, tgt_type = :Form1, op_names = [:⋆, :⋆₁⁻¹]),

  # Rules for Δ
  (src_type = :Form0, tgt_type = :Form0, op_names = [:Δ, :Δ₀]),
  (src_type = :Form1, tgt_type = :Form1, op_names = [:Δ, :Δ₁]),

  # Rules for δ
  (src_type = :Form1, tgt_type = :Form0, op_names = [:δ, :δ₁]),

  # Rules for negation
  (src_type = :Form0, tgt_type = :Form0, op_names = [:neg, :(-)]),
  (src_type = :Form1, tgt_type = :Form1, op_names = [:neg, :(-)])]

op2_inf_rules_1D = [
  # Rules for ∧₀₀, ∧₁₀, ∧₀₁
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :Form0, op_names = [:∧, :∧₀₀]),
  (proj1_type = :Form1, proj2_type = :Form0, res_type = :Form1, op_names = [:∧, :∧₁₀]),
  (proj1_type = :Form0, proj2_type = :Form1, res_type = :Form1, op_names = [:∧, :∧₀₁]),

  # Rules for L₀, L₁
  (proj1_type = :Form1, proj2_type = :Form0, res_type = :Form0, op_names = [:L, :L₀]),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form1, op_names = [:L, :L₁]),    

  # Rules for i₁
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form0, op_names = [:i, :i₁]),

  # Rules for divison and multiplication
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :Form0, op_names = [:./, :.*]),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form1, op_names = [:./, :.*]),
  (proj1_type = :Form2, proj2_type = :Form2, res_type = :Form2, op_names = [:./, :.*]),

  # This parameter type inference might be wrong, depending on what the user gives as a parameter
  (proj1_type = :Parameter, proj2_type = :Form0, res_type = :Form0, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Parameter, proj2_type = :Form1, res_type = :Form1, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Parameter, proj2_type = :Form2, res_type = :Form2, op_names = [:/, :./, :*, :.*]),

  (proj1_type = :Form0, proj2_type = :Parameter, res_type = :Form0, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Form1, proj2_type = :Parameter, res_type = :Form1, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Form2, proj2_type = :Parameter, res_type = :Form2, op_names = [:/, :./, :*, :.*]),
  
  (proj1_type = :Literal, proj2_type = :Form0, res_type = :Form0, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Literal, proj2_type = :Form1, res_type = :Form1, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Literal, proj2_type = :Form2, res_type = :Form2, op_names = [:/, :./, :*, :.*]),
  
  (proj1_type = :Form0, proj2_type = :Literal, res_type = :Form0, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Form1, proj2_type = :Literal, res_type = :Form1, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Form2, proj2_type = :Literal, res_type = :Form2, op_names = [:/, :./, :*, :.*])]

  
"""
These are the default rules used to do type inference in the 2D exterior calculus.
"""
op1_inf_rules_2D = [
  # Rules for ∂ₜ
  (src_type = :Form0, tgt_type = :Form0, op_names = [:∂ₜ]),
  (src_type = :Form1, tgt_type = :Form1, op_names = [:∂ₜ]),
  (src_type = :Form2, tgt_type = :Form2, op_names = [:∂ₜ]),

  # Rules for d
  (src_type = :Form0, tgt_type = :Form1, op_names = [:d, :d₀]),
  (src_type = :Form1, tgt_type = :Form2, op_names = [:d, :d₁]),
  (src_type = :DualForm0, tgt_type = :DualForm1, op_names = [:d, :dual_d₀, :d̃₀]),
  (src_type = :DualForm1, tgt_type = :DualForm2, op_names = [:d, :dual_d₁, :d̃₁]),

  # Rules for ⋆
  (src_type = :Form0, tgt_type = :DualForm2, op_names = [:⋆, :⋆₀]),
  (src_type = :Form1, tgt_type = :DualForm1, op_names = [:⋆, :⋆₁]),
  (src_type = :Form2, tgt_type = :DualForm0, op_names = [:⋆, :⋆₂]),

  (src_type = :DualForm2, tgt_type = :Form0, op_names = [:⋆, :⋆₀⁻¹]),
  (src_type = :DualForm1, tgt_type = :Form1, op_names = [:⋆, :⋆₁⁻¹]),
  (src_type = :DualForm0, tgt_type = :Form2, op_names = [:⋆, :⋆₂⁻¹]),

  # Rules for Δ
  (src_type = :Form0, tgt_type = :Form0, op_names = [:Δ, :Δ₀]),
  (src_type = :Form1, tgt_type = :Form1, op_names = [:Δ, :Δ₁]),
  (src_type = :Form2, tgt_type = :Form2, op_names = [:Δ, :Δ₂]),

  # Rules for δ
  (src_type = :Form1, tgt_type = :Form0, op_names = [:δ, :δ₁]),
  (src_type = :Form2, tgt_type = :Form1, op_names = [:δ, :δ₂]),

  # Rules for negation
  (src_type = :Form0, tgt_type = :Form0, op_names = [:neg, :(-)]),
  (src_type = :Form1, tgt_type = :Form1, op_names = [:neg, :(-)]),
  (src_type = :Form2, tgt_type = :Form2, op_names = [:neg, :(-)])]
    
op2_inf_rules_2D = vcat(op2_inf_rules_1D, [
  # Rules for ∧₁₁, ∧₂₀, ∧₀₂
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form2, op_names = [:∧, :∧₁₁]),
  (proj1_type = :Form2, proj2_type = :Form0, res_type = :Form2, op_names = [:∧, :∧₂₀]),
  (proj1_type = :Form0, proj2_type = :Form2, res_type = :Form2, op_names = [:∧, :∧₀₂]),

  # Rules for L₂
  (proj1_type = :Form1, proj2_type = :Form2, res_type = :Form2, op_names = [:L, :L₂]),    

  # Rules for i₁
  (proj1_type = :Form1, proj2_type = :Form2, res_type = :Form1, op_names = [:i, :i₁])])
  
function apply_inference_rule_op1!(d::SummationDecapode, op1_id, rule)
  type_src = d[d[op1_id, :src], :type]
  type_tgt = d[d[op1_id, :tgt], :type]

  if(type_src != :infer && type_tgt != :infer)
    return false
  end

  score_src = (rule.src_type == type_src)
  score_tgt = (rule.tgt_type == type_tgt)
  check_op = (d[op1_id, :op1] in rule.op_names)

  if(check_op && (score_src + score_tgt == 1))
    d[d[op1_id, :src], :type] = rule.src_type
    d[d[op1_id, :tgt], :type] = rule.tgt_type
    return true
  end

  return false
end

function apply_inference_rule_op2!(d::SummationDecapode, op2_id, rule)
  type_proj1 = d[d[op2_id, :proj1], :type]
  type_proj2 = d[d[op2_id, :proj2], :type]
  type_res = d[d[op2_id, :res], :type]

  if(type_proj1 != :infer && type_proj2 != :infer && type_res != :infer)
    return false
  end

  score_proj1 = (rule.proj1_type == type_proj1)
  score_proj2 = (rule.proj2_type == type_proj2)
  score_res = (rule.res_type == type_res)
  check_op = (d[op2_id, :op2] in rule.op_names)

  if(check_op && (score_proj1 + score_proj2 + score_res == 2))
    d[d[op2_id, :proj1], :type] = rule.proj1_type
    d[d[op2_id, :proj2], :type] = rule.proj2_type
    d[d[op2_id, :res], :type] = rule.res_type
    return true
  end

  return false
end

"""
  function infer_summands_and_summations!(d::SummationDecapode)

"""
function infer_summands_and_summations!(d::SummationDecapode)
  # Note that we are not doing any type checking here!
  # i.e. We are not checking for this: [Form0, Form1, Form0].
  applied = false
  for Σ_idx in parts(d, :Σ)
    summands = d[:summand][incident(d, Σ_idx, :summation)]
    sum = d[:sum][Σ_idx]
    idxs = [summands; sum]
    types = d[:type][idxs]
    all(t != :infer for t in types) && continue # We need not infer
    all(t == :infer for t in types) && continue # We can  not infer
    inferred_type = types[findfirst(!=(:infer), types)]
    to_infer_idxs = filter(i -> d[:type][i] == :infer, idxs)
    d[:type][to_infer_idxs] .= inferred_type
    applied = true
  end
  return applied
end

# TODO: Although the big-O complexity is the same, it might be more efficent on
# average to iterate over edges then rules, instead of rules then edges. This
# might result in more un-maintainable code. If you implement this, you might
# also want to make the rules keys in a Dict.
# It also might be more efficient on average to instead iterate over variables.
"""
  function infer_types!(d::SummationDecapode, op1_rules::Vector{NamedTuple{(:src_type, :tgt_type, :replacement_type, :op), NTuple{4, Symbol}}})

Infer types of Vars given rules wherein one type is known and the other not.
"""
function infer_types!(d::SummationDecapode, op1_rules::Vector{NamedTuple{(:src_type, :tgt_type, :op_names), Tuple{Symbol, Symbol, Vector{Symbol}}}}, op2_rules::Vector{NamedTuple{(:proj1_type, :proj2_type, :res_type, :op_names), Tuple{Symbol, Symbol, Symbol, Vector{Symbol}}}})
  # This is an optimization so we do not "visit" a row which has no infer types.
  # It could be deleted if found to be not worth maintainability tradeoff.
  types_known_op1 = ones(Bool, nparts(d, :Op1))
  types_known_op1[incident(d, :infer, [:src, :type])] .= false
  types_known_op1[incident(d, :infer, [:tgt, :type])] .= false

  #types_known_op1 = zeros(Bool, nparts(d, :Op1))
  #types_known_op1[incident(d, :infer, [:src, :type])] .= false
  #types_known_op1[incident(d, :infer, [:tgt, :type])] .= false

  types_known_op2 = zeros(Bool, nparts(d, :Op2))
  types_known_op2[incident(d, :infer, [:proj1, :type])] .= false
  types_known_op2[incident(d, :infer, [:proj2, :type])] .= false
  types_known_op2[incident(d, :infer, [:res, :type])] .= false

  while true
    applied = false
    
    for rule in op1_rules
      for op1_idx in parts(d, :Op1)
        types_known_op1[op1_idx] && continue

        this_applied = apply_inference_rule_op1!(d, op1_idx, rule)

        types_known_op1[op1_idx] = this_applied
        applied = applied || this_applied
      end
    end

    for rule in op2_rules
      for op2_idx in parts(d, :Op2)
        types_known_op2[op2_idx] && continue

        this_applied = apply_inference_rule_op2!(d, op2_idx, rule)

        types_known_op2[op2_idx] = this_applied
        applied = applied || this_applied
      end
    end

    applied = applied || infer_summands_and_summations!(d)
    applied || break # Break if no rules were applied.
  end 

  d
end

# TODO: When SummationDecapodes are annotated with the degree of their space,
# use dispatch to choose the correct set of rules.
infer_types!(d::SummationDecapode) =
  infer_types!(d, op1_inf_rules_2D, op2_inf_rules_2D)

# TODO: You could write a method which auto-generates these rules given degree N.

"""
These are the default rules used to do function resolution in the 1D exterior calculus.
"""
op1_res_rules_1D = [
  # Rules for d.
  (src_type = :Form0, tgt_type = :Form1, resolved_name = :d₀, op = :d),
  (src_type = :DualForm0, tgt_type = :DualForm1, resolved_name = :dual_d₀, op = :d),
  # Rules for ⋆.
  (src_type = :Form0, tgt_type = :DualForm1, resolved_name = :⋆₀, op = :⋆),
  (src_type = :Form1, tgt_type = :DualForm0, resolved_name = :⋆₁, op = :⋆),
  (src_type = :DualForm1, tgt_type = :Form0, resolved_name = :⋆₀⁻¹, op = :⋆),
  (src_type = :DualForm0, tgt_type = :Form1, resolved_name = :⋆₁⁻¹, op = :⋆),
  # Rules for δ.
  (src_type = :Form1, tgt_type = :Form0, resolved_name = :δ₁, op = :δ)]

# We merge 1D and 2D rules since it seems op2 rules are metric-free. If
# this assumption is false, this needs to change.
op2_res_rules_1D = [
  # Rules for ∧.
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :Form0, resolved_name = :∧₀₀, op = :∧),
  (proj1_type = :Form1, proj2_type = :Form0, res_type = :Form1, resolved_name = :∧₁₀, op = :∧),
  (proj1_type = :Form0, proj2_type = :Form1, res_type = :Form1, resolved_name = :∧₀₁, op = :∧),
  # Rules for L.
  (proj1_type = :Form1, proj2_type = :Form0, res_type = :Form0, resolved_name = :L₀, op = :L),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form1, resolved_name = :L₁, op = :L),
  # Rules for i.
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form0, resolved_name = :i₁, op = :i)]


"""
These are the default rules used to do function resolution in the 2D exterior calculus.
"""
op1_res_rules_2D = [
  # Rules for d.
  (src_type = :Form0, tgt_type = :Form1, resolved_name = :d₀, op = :d),
  (src_type = :Form1, tgt_type = :Form2, resolved_name = :d₁, op = :d),
  (src_type = :DualForm0, tgt_type = :DualForm1, resolved_name = :dual_d₀, op = :d),
  (src_type = :DualForm1, tgt_type = :DualForm2, resolved_name = :dual_d₁, op = :d),
  # Rules for ⋆.
  (src_type = :Form0, tgt_type = :DualForm2, resolved_name = :⋆₀, op = :⋆),
  (src_type = :Form1, tgt_type = :DualForm1, resolved_name = :⋆₁, op = :⋆),
  (src_type = :Form2, tgt_type = :DualForm0, resolved_name = :⋆₂, op = :⋆),
  (src_type = :DualForm2, tgt_type = :Form0, resolved_name = :⋆₀⁻¹, op = :⋆),
  (src_type = :DualForm1, tgt_type = :Form1, resolved_name = :⋆₁⁻¹, op = :⋆),
  (src_type = :DualForm0, tgt_type = :Form2, resolved_name = :⋆₂⁻¹, op = :⋆),
  # Rules for δ.
  (src_type = :Form2, tgt_type = :Form1, resolved_name = :δ₂, op = :δ),
  (src_type = :Form1, tgt_type = :Form0, resolved_name = :δ₁, op = :δ),
  # Rules for ∇².
  (src_type = :Form0, tgt_type = :Form0, resolved_name = :∇²₀, op = :∇²),
  (src_type = :Form1, tgt_type = :Form1, resolved_name = :∇²₁, op = :∇²),
  (src_type = :Form2, tgt_type = :Form2, resolved_name = :∇²₂, op = :∇²),
  # Rules for Δ.
  (src_type = :Form0, tgt_type = :Form0, resolved_name = :Δ₀, op = :Δ),
  (src_type = :Form1, tgt_type = :Form1, resolved_name = :Δ₁, op = :Δ),
  (src_type = :Form1, tgt_type = :Form1, resolved_name = :Δ₂, op = :Δ)]

# We merge 1D and 2D rules directly here since it seems op2 rules
# are metric-free. If this assumption is false, this needs to change.
op2_res_rules_2D = vcat(op2_res_rules_1D, [
  # Rules for ∧.
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form2, resolved_name = :∧₁₁, op = :∧),
  (proj1_type = :Form2, proj2_type = :Form0, res_type = :Form2, resolved_name = :∧₂₀, op = :∧),
  (proj1_type = :Form0, proj2_type = :Form2, res_type = :Form2, resolved_name = :∧₀₂, op = :∧),
  # Rules for L.
  (proj1_type = :Form1, proj2_type = :Form2, res_type = :Form2, resolved_name = :L₂, op = :L),
  # Rules for i.
  (proj1_type = :Form1, proj2_type = :Form2, res_type = :Form1, resolved_name = :i₂, op = :i)])
  
"""
  function resolve_overloads!(d::SummationDecapode, op1_rules::Vector{NamedTuple{(:src_type, :tgt_type, :resolved_name, :op), NTuple{4, Symbol}}})

Resolve function overloads based on types of src and tgt.
"""
function resolve_overloads!(d::SummationDecapode, op1_rules::Vector{NamedTuple{(:src_type, :tgt_type, :resolved_name, :op), NTuple{4, Symbol}}}, op2_rules::Vector{NamedTuple{(:proj1_type, :proj2_type, :res_type, :resolved_name, :op), NTuple{5, Symbol}}})
  for op1_idx in parts(d, :Op1)
    src = d[:src][op1_idx]; tgt = d[:tgt][op1_idx]; op1 = d[:op1][op1_idx]
    src_type = d[:type][src]; tgt_type = d[:type][tgt]
    for rule in op1_rules
      if op1 == rule[:op] && src_type == rule[:src_type] && tgt_type == rule[:tgt_type]
        d[:op1][op1_idx] = rule[:resolved_name]
        break
      end
    end
  end

  for op2_idx in parts(d, :Op2)
    proj1 = d[:proj1][op2_idx]; proj2 = d[:proj2][op2_idx]; res = d[:res][op2_idx]; op2 = d[:op2][op2_idx]
    proj1_type = d[:type][proj1]; proj2_type = d[:type][proj2]; res_type = d[:type][res]
    for rule in op2_rules
      if op2 == rule[:op] && proj1_type == rule[:proj1_type] && proj2_type == rule[:proj2_type] && res_type == rule[:res_type]
        d[:op2][op2_idx] = rule[:resolved_name]
        break
      end
    end
  end

  d
end

# TODO: When SummationDecapodes are annotated with the degree of their space,
# use dispatch to choose the correct set of rules.
resolve_overloads!(d::SummationDecapode) =
  resolve_overloads!(d, op1_res_rules_2D, op2_res_rules_2D)

