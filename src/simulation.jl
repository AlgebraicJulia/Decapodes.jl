using Base.Iterators
using CombinatorialSpaces
using ComponentArrays
using LinearAlgebra
using MLStyle
using PreallocationTools

const gensim_in_place_stub = Symbol("GenSim-M")

abstract type GenerationTarget end

struct CPUTarget <: GenerationTarget end
struct CUDATarget <: GenerationTarget end

abstract type AbstractCall end

struct UnaryCall <: AbstractCall
  operator
  equality
  input
  output
end

# TODO: Add back support for contract operators
Base.Expr(c::UnaryCall) = begin
  operator = c.operator
  #= if isa(operator, AbstractArray)
    operator = Expr(:call, :∘, reverse(operator)...)
  end =#
  if(c.equality == :.=)
    if(operator == add_inplace_stub(:⋆₁⁻¹)) # Since inverse hodge Geo is a solver
      Expr(:call, c.operator, c.output, c.input)
    elseif(operator == :.-)
      Expr(c.equality, c.output, Expr(:call, operator, c.input))
    else
      Expr(:call, :mul!, c.output, operator, c.input)
    end
  else
    Expr(c.equality, c.output, Expr(:call, operator, c.input))
  end
end

struct BinaryCall <: AbstractCall
  operator
  equality
  input1
  input2
  output
end

# TODO: After getting rid of AppCirc2, do we need this check?
Base.Expr(c::BinaryCall) = begin
  #= if isa(c.operator, AbstractArray)
    operator = :(compose($(c.operator)))
  end =#

  # These operators can be done in-place
  if(c.equality == :.= && get_stub(c.operator) == gensim_in_place_stub)
    return Expr(:call, c.operator, c.output, c.input1, c.input2)
  end
  return Expr(c.equality, c.output, Expr(:call, c.operator, c.input1, c.input2))
end

struct VarargsCall <: AbstractCall
  operator
  equality
  inputs
  output
end

Base.Expr(c::VarargsCall) = begin
  #= if isa(c.operator, AbstractArray)
    operator = :(compose($(c.operator)))
  end =#
  return Expr(c.equality, c.output, Expr(:call, c.operator, c.inputs...))
end

struct AllocVecCall <: AbstractCall
  name
  form
  dimension
  T
  code_target
end

struct AllocVecCallError <: Exception
  c::AllocVecCall
end

# TODO: There are likely better ways of dispatching on dimension instead of
# storing it inside an AllocVecCall.
Base.Expr(c::AllocVecCall) = begin
  resolved_form = @match (c.form, c.dimension) begin
    (:Form0, 2) => :V
    (:Form1, 2) => :E
    (:Form2, 2) => :Tri
    (:DualForm0, 2) => :Tri
    (:DualForm1, 2) => :E
    (:DualForm2, 2) => :V

    (:Form0, 1) => :V
    (:Form1, 1) => :E
    (:DualForm0, 1) => :E
    (:DualForm1, 1) => :V
    _ => throw(AllocVecCallError(c))
  end

  hook_AVC_caching(c, resolved_form, c.code_target)
end

function hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::CPUTarget)
  :($(Symbol(:__,c.name)) = Decapodes.FixedSizeDiffCache(Vector{$(c.T)}(undef, nparts(mesh, $(QuoteNode(resolved_form))))))
end

# TODO: Allow user to overload these hooks with user-defined code_target
function hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::CUDATarget)
  :($(c.name) = CuVector{$(c.T)}(undef, nparts(mesh, $(QuoteNode(resolved_form)))))
end

#= function get_form_number(d::SummationDecapode, var_id::Int)
  type = d[var_id, :type]
  if(type == :Form0)
    return 0
  elseif(type == :Form1)
    return 1
  elseif(type == :Form2)
    return 2
  end
  return -1
end

# WARNING: This may not work if names are not unique, use above instead
function get_form_number(d::SummationDecapode, var_name::Symbol)
var_id = first(incident(d, var_name, :name))
return get_form_number(d, var_id)
end =#

function is_form(d::SummationDecapode, var_id::Int)
  type = d[var_id, :type]
  return (type == :Form0 || type == :Form1 || type == :Form2 ||
    type == :DualForm0 || type == :DualForm1 || type == :DualForm2)
end

is_form(d::SummationDecapode, var_name::Symbol) = is_form(d, first(incident(d, var_name, :name)))

is_literal(d::SummationDecapode, var_id::Int) = (d[var_id, :type] == :Literal)
is_literal(d::SummationDecapode, var_name::Symbol) = is_literal(d, first(incident(d, var_name, :name)))

is_infer(d::SummationDecapode, var_id::Int) = (d[var_id, :type] == :infer)
is_infer(d::SummationDecapode, var_name::Symbol) = is_infer(d, first(incident(d, var_name, :name)))

add_stub(stub_name::Symbol, var_name::Symbol) = return Symbol("$(stub_name)_$(var_name)")

function get_stub(var_name::Symbol)
  var_str = String(var_name)
  idx = findfirst("_", var_str)
  if(isnothing(idx) || first(idx) == 1)
    return nothing
  end
  return Symbol(var_str[begin:first(idx) - 1])
end

add_inplace_stub(var_name::Symbol) = add_stub(gensim_in_place_stub, var_name)

const ARITHMETIC_OPS = Set([:+, :*, :-, :/, :.+, :.*, :.-, :./, :^, :.^, :.>, :.<, :.≤, :.≥])

"""
    compile_env(d::AbstractNamedDecapode, dec_matrices::Vector{Symbol}, con_dec_operators::Set{Symbol}, code_target::GenerationTarget = CPUTarget())

This creates the symbol to function linking for the simulation output. Those run through the `default_dec` backend
expect both an in-place and an out-of-place variant in that order. User defined operations only support out-of-place.
"""
function compile_env(d::AbstractNamedDecapode, dec_matrices::Vector{Symbol}, con_dec_operators::Set{Symbol}, code_target::GenerationTarget = CPUTarget())
  defined_ops = deepcopy(con_dec_operators)

  defs = quote end

  for op in dec_matrices
    if op in defined_ops
      continue
    end

    quote_op = QuoteNode(op)
    mat_op = add_stub(gensim_in_place_stub, op)

    # TODO: Add support for user-defined code targets
    default_generation = @match code_target begin
      ::CPUTarget => :default_dec_matrix_generate
      ::CUDATarget => :default_dec_cu_matrix_generate
      _ => error("Provided code target $(code_target) is not yet supported in simulations")
    end

    def = :(($mat_op, $op) = $(default_generation)(mesh, $quote_op, hodge))
    push!(defs.args, def)

    push!(defined_ops, op)
  end

  # Add in user-defined operations
  for op in vcat(d[:op1], d[:op2])
    if op == DerivOp || op in defined_ops || op in ARITHMETIC_OPS
      continue
    end
    ops = QuoteNode(op)
    def = :($op = operators(mesh, $ops))
    push!(defs.args, def)

    push!(defined_ops, op)
  end

  return defs
end

function compile_var(alloc_vectors::Vector{AllocVecCall})
  return quote $(Expr.(alloc_vectors)...) end
end

# This is the block of parameter setting inside f
# TODO: Pass this an extra type parameter that sets the size of the Floats
get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol}) = get_vars_code(d, vars, Float64)

"""
    get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol}, ::Type{stateeltype}) where stateeltype

This initalizes all input variables according to their Decapodes type.
"""
function get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol}, ::Type{stateeltype}) where stateeltype
  stmts = map(vars) do s
    ssymbl = QuoteNode(s)
    if all(d[incident(d, s, :name) , :type] .== :Constant)
      :($s = p.$s)
    elseif all(d[incident(d, s, :name) , :type] .== :Parameter)
      :($s = (p.$s)(t))
    elseif all(d[incident(d, s, :name) , :type] .== :Literal)
      # Literals don't need assignments, because they are literals, but we stored them as Symbols.
      # #TODO: we should fix that upstream so that we don't need this.
      :($s = $(parse(stateeltype, String(s))))
    else
      # TODO: If names are not unique, then the type is assumed to be a
      # form for all of the vars sharing a same name.
      :($s = u.$s)
    end
  end
  return quote $(stmts...) end
end

"""
    set_tanvars_code(d::AbstractNamedDecapode)

This is
"""
function set_tanvars_code(d::AbstractNamedDecapode)
  tanvars = [(d[e, [:src,:name]], d[e, [:tgt,:name]]) for e in incident(d, :∂ₜ, :op1)]
  stmts = map(tanvars) do (s,t)
    ssymb = QuoteNode(s)
    :(getproperty(du, $ssymb) .= $t)
  end
  return stmts
end

function compile(d::SummationDecapode, inputs::Vector, alloc_vectors::Vector{AllocVecCall}, optimizable_dec_operators::Set{Symbol}; dimension=2, stateeltype=Float64, code_target=CPUTarget())
  # Get the Vars of the inputs (probably state Vars).
  visited_Var = falses(nparts(d, :Var))

  # input_numbers = incident(d, inputs, :name)
  input_numbers = reduce(vcat, incident(d, inputs, :name))
  # visited_Var[collect(flatten(input_numbers))] .= true
  visited_Var[input_numbers] .= true
  visited_Var[incident(d, :Literal, :type)] .= true

  visited_1 = falses(nparts(d, :Op1))
  visited_2 = falses(nparts(d, :Op2))
  visited_Σ = falses(nparts(d, :Σ))

  promote_arithmetic_map = Dict(:(+) => :.+, :(-) => :.-, :(*) => :.*, :(/) => :./, :(^) => :.^, :(=) => :.=,
                                :.+ => :.+, :.- => :.-, :.* => :.*, :./ => :./, :.^ => :.^, :.= => :.=)

  # FIXME: this is a quadratic implementation of topological_sort inlined in here.
  op_order = []
  for _ in 1:(nparts(d, :Op1) + nparts(d,:Op2) + nparts(d, :Σ))
    for op in parts(d, :Op1)
      s = d[op, :src]
      if !visited_1[op] && visited_Var[s]
        # skip the derivative edges
        operator = d[op, :op1]
        t = d[op, :tgt]
        visited_1[op] = true
        if operator == DerivOp
          continue
        end

        equality = :(=)

        sname = d[s, :name]
        tname = d[t, :name]

        # TODO: Check to see if this is a DEC operator
        if(operator in optimizable_dec_operators)
          # push!(dec_matrices, operator)
          if(is_form(d, t))
            equality = promote_arithmetic_map[equality]
            operator = add_stub(gensim_in_place_stub, operator)

            push!(alloc_vectors, AllocVecCall(tname, d[t, :type], dimension, stateeltype, code_target))
          end
        elseif(operator == :(-) || operator == :.-)
          if(is_form(d, t))
            equality = promote_arithmetic_map[equality]
            operator = promote_arithmetic_map[operator]

            push!(alloc_vectors, AllocVecCall(tname, d[t, :type], dimension, stateeltype, code_target))
          end
        end

        visited_Var[t] = true
        c = UnaryCall(operator, equality, sname, tname)
        push!(op_order, c)
      end
    end

    for op in parts(d, :Op2)
      arg1 = d[op, :proj1]
      arg2 = d[op, :proj2]
      if !visited_2[op] && visited_Var[arg1] && visited_Var[arg2]
        r = d[op, :res]
        a1name = d[arg1, :name]
        a2name = d[arg2, :name]
        rname  = d[r, :name]

        operator = d[op, :op2]
        equality = :(=)

        # This is meant for wedge products
        #= if(operator in optimizable_dec_operators)
          push!(dec_matrices, operator)
        end =#

        # TODO: Check to make sure that this logic never breaks
        if(is_form(d, r))
          if(operator == :(+) || operator == :(-) || operator == :.+ || operator == :.-)
            operator = promote_arithmetic_map[operator]
            equality = promote_arithmetic_map[equality]
            push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))

          # TODO: Do we want to support the ability of a user to use the backslash operator?
          elseif(operator == :(*) || operator == :(/) || operator == :.* || operator == :./)
            # WARNING: This part may break if we add more compiler types that have different
            # operations for basic and broadcast modes, e.g. matrix multiplication vs broadcast
            if(!is_infer(d, arg1) && !is_infer(d, arg2))
              operator = promote_arithmetic_map[operator]
              equality = promote_arithmetic_map[equality]
              push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))
            end
          elseif(operator in optimizable_dec_operators)
            operator = add_stub(gensim_in_place_stub, operator)
            equality = promote_arithmetic_map[equality]
            push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))
          end
        end

        # TODO: Clean this in another PR (with a @match maybe).
        if(operator == :(*))
          operator = promote_arithmetic_map[operator]
        end
        if(operator == :(-))
          operator = promote_arithmetic_map[operator]
        end
        if(operator == :(/))
          operator = promote_arithmetic_map[operator]
        end
        if(operator == :(^))
          operator = promote_arithmetic_map[operator]
        end

        visited_2[op] = true
        visited_Var[r] = true
        c = BinaryCall(operator, equality, a1name, a2name, rname)
        push!(op_order, c)
      end
    end

    for op in parts(d, :Σ)
      args = subpart(d, incident(d, op, :summation), :summand)
      if !visited_Σ[op] && all(visited_Var[args])
        r = d[op, :sum]
        argnames = d[args, :name]
        rname  = d[r, :name]

        operator = :(+)
        equality = :(=)

        # If result is a known form, broadcast addition
        # TODO: Also need to tell handler to prealloc
        if(is_form(d, r))
          operator = promote_arithmetic_map[operator]
          equality = promote_arithmetic_map[equality]
          push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))
        end

        operator = :(.+)

        visited_Σ[op] = true
        visited_Var[r] = true
        c = VarargsCall(operator, equality, argnames, rname)
        push!(op_order, c)
      end
    end
  end

  cache_exprs = []
  if(code_target isa CPUTarget)
    cache_exprs = map(alloc_vectors) do vec
      :($(vec.name) = (Decapodes.get_tmp($(Symbol(:__,vec.name)), u)))
    end
  end

  eq_exprs = map(Expr, op_order)
  vcat(cache_exprs,eq_exprs)
end

# TODO: Add more specific types later for optimization
function resolve_types_compiler!(d::SummationDecapode)
  d[:type] = map(d[:type]) do x
    if(x == :Constant || x == :Parameter)
      return :infer
    end
    return x
  end
end

function replace_names_compiler!(d::SummationDecapode)
  dec_op1 = Pair{Symbol, Any}[]
  dec_op2 = Pair{Symbol, Symbol}[(:∧₀₀ => :.*)]
  replace_names!(d, dec_op1, dec_op2)
end

function infer_overload_compiler!(d::SummationDecapode, dimension::Int)
  if(dimension == 1)
    infer_types!(d, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(d, op1_res_rules_1D, op2_res_rules_1D)
  elseif(dimension == 2)
    infer_types!(d, op1_inf_rules_2D, op2_inf_rules_2D)
    resolve_overloads!(d, op1_res_rules_2D, op2_res_rules_2D)
  end
end

function init_dec_matrices!(d::SummationDecapode, dec_matrices::Vector{Symbol}, optimizable_dec_operators::Set{Symbol})

  for op1_name in d[:op1]
    if(op1_name ∈ optimizable_dec_operators)
      push!(dec_matrices, op1_name)
    end
  end

  for op2_name in d[:op2]
    if(op2_name ∈ optimizable_dec_operators)
      push!(dec_matrices, op2_name)
    end
  end
end

function link_contract_operators(d::SummationDecapode, con_dec_operators::Set{Symbol}, code_target::GenerationTarget, stateeltype::DataType)

  contract_defs = quote end

  compute_to_name = Dict()
  curr_id = 1

  for op1_id in parts(d, :Op1)
    op1_name = d[op1_id, :op1]
    if isa(op1_name, AbstractArray)
      computation = reverse!(map(x -> add_inplace_stub(x), op1_name))
      compute_key = join(computation, " * ")

      computation_name = get(compute_to_name, compute_key, :Error)
      if(computation_name == :Error)
        computation_name = add_stub(Symbol("GenSim-ConMat"), Symbol(curr_id))
        get!(compute_to_name, compute_key, computation_name)
        push!(con_dec_operators, computation_name)

        expr_line = hook_LCO_inplace(computation_name, computation, code_target, stateeltype)
        push!(contract_defs.args, expr_line)

        expr_line = Expr(Symbol("="), computation_name, Expr(Symbol("->"), :x, Expr(:call, :*, add_inplace_stub(computation_name), :x)))
        push!(contract_defs.args, expr_line)

        curr_id += 1
      end

      d[op1_id, :op1] = computation_name
    end
  end

  contract_defs
end

# TODO: Allow user to overload these hooks with user-defined code_target
function hook_LCO_inplace(computation_name, computation, ::CPUTarget, float_type::DataType)
  return :($(add_inplace_stub(computation_name)) = $(Expr(:call, :*, computation...)))
end

# Turn multiplication right-associative to prevent CUDA.jl error
# TODO: Can turn this into left-associative
function generate_parentheses_multipy(list)
  if(length(list) == 1)
      return list[1]
  else
      return Expr(:call, :*, list[1], generate_parentheses_multipy(list[2:end]))
  end
end

# Adapt this to also write Diagonal Matrices as CuVectors, sparsifying diagonal matrices slows down computations 
function hook_LCO_inplace(computation_name, computation, ::CUDATarget, float_type::DataType)
  # return :($(add_inplace_stub(computation_name)) = CUDA.CUSPARSE.CuSparseMatrixCSC{$(float_type)}($(Expr(:call, :*, computation...))))
  return :($(add_inplace_stub(computation_name)) = $(generate_parentheses_multipy(computation)))
end

function gensim(user_d::AbstractNamedDecapode, input_vars; dimension::Int=2, stateeltype = Float64, code_target = CPUTarget())
  recognize_types(user_d)

  # Makes copy
  d′ = expand_operators(user_d)

  dec_matrices = Vector{Symbol}();
  alloc_vectors = Vector{AllocVecCall}();

  vars = get_vars_code(d′, input_vars, stateeltype)
  tars = set_tanvars_code(d′)

  # We need to run this after we grab the constants and parameters out
  infer_overload_compiler!(d′, dimension)
  resolve_types_compiler!(d′)
  infer_overload_compiler!(d′, dimension)

  # This should probably be followed by an expand_operators
  replace_names_compiler!(d′)
  open_operators!(d′, dimension = dimension)
  infer_overload_compiler!(d′, dimension)

  # This will generate all of the fundemental DEC operators present
  optimizable_dec_operators = Set([:⋆₀, :⋆₁, :⋆₂, :⋆₀⁻¹, :⋆₂⁻¹,
                                  :d₀, :d₁, :dual_d₀, :d̃₀, :dual_d₁, :d̃₁])
  extra_dec_operators = Set([:⋆₁⁻¹, :∧₀₁, :∧₁₀, :∧₁₁, :∧₀₂, :∧₂₀])

  init_dec_matrices!(d′, dec_matrices, union(optimizable_dec_operators, extra_dec_operators))

  # This contracts matrices together into a single matrix
  contracted_dec_operators = Set{Symbol}();
  contract_operators!(d′, allowable_ops = optimizable_dec_operators)
  cont_defs = link_contract_operators(d′, contracted_dec_operators, code_target, stateeltype)

  union!(optimizable_dec_operators, contracted_dec_operators, extra_dec_operators)

  # Compilation of the simulation
  equations = compile(d′, input_vars, alloc_vectors, optimizable_dec_operators, dimension=dimension, stateeltype=stateeltype, code_target=code_target)

  func_defs = compile_env(d′, dec_matrices, contracted_dec_operators, code_target)
  vect_defs = compile_var(alloc_vectors)

  quote
    (mesh, operators, hodge=GeometricHodge()) -> begin
      $func_defs
      $cont_defs
      $vect_defs
      f(du, u, p, t) = begin
        $vars
        $(equations...)
        $(tars...)
        return nothing
      end;
    end
  end
end

gensim(c::Collage; dimension::Int=2) =
gensim(collate(c); dimension=dimension)

"""    function gensim(d::AbstractNamedDecapode; dimension::Int=2)

Generate a simulation function from the given Decapode. The returned function can then be combined with a mesh and a function describing function mappings to return a simulator to be passed to `solve`.
"""
gensim(d::AbstractNamedDecapode; dimension::Int=2, stateeltype = Float64, code_target = CPUTarget()) =
  gensim(d, vcat(collect(infer_state_names(d)), d[incident(d, :Literal, :type), :name]), dimension=dimension, stateeltype=stateeltype, code_target=code_target)

evalsim(d::AbstractNamedDecapode; dimension::Int=2, stateeltype = Float64, code_target = CPUTarget()) =
  eval(gensim(d, dimension=dimension, stateeltype=stateeltype, code_target=code_target))
evalsim(d::AbstractNamedDecapode, input_vars; dimension::Int=2, stateeltype = Float64, code_target = CPUTarget()) =
  eval(gensim(d, input_vars, dimension=dimension, stateeltype=stateeltype, code_target=code_target))

"""
function find_unreachable_tvars(d)

Determine if the given Decapode can be compiled to an explicit time-stepping
simulation. Use the simple check that one can traverse the Decapode starting
from the state variables and reach all TVars (ignoring ∂ₜ edges).
"""
function find_unreachable_tvars(d)
  states = infer_states(d)
  visited = falses(nparts(d, :Var))
  s = Stack{Int64}()
  foreach(v -> push!(s, v), states)
  while true
    curr_var = pop!(s)
    visited[curr_var] = true
    op1s = incident(d, curr_var, :src)
    filter!(x -> d[x, :op1] == :∂ₜ, op1s)
    for op1 in op1s
      (!visited[d[op1, :tgt]]) && push!(s, d[op1, :tgt])
    end
  end
  TVars = d[:incl]
end
