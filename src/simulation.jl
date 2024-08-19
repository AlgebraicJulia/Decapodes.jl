using CombinatorialSpaces
using ComponentArrays
using LinearAlgebra
using MLStyle
using PreallocationTools

const GENSIM_INPLACE_STUB = Symbol("GenSim-M")
const NO_STUB_RETURN = Symbol("NOSTUB")

abstract type AbstractGenerationTarget end

abstract type CPUBackend <: AbstractGenerationTarget end
abstract type CUDABackend <: AbstractGenerationTarget end

struct CPUTarget <: CPUBackend end
struct CUDATarget <: CUDABackend end

# TODO: Make it so AbstractCall code terminates into an error, not into default code
abstract type AbstractCall end

struct InvalidCallException <: Exception end

Base.showerror(io::IO, e::InvalidCallException) = print(io, "Compiler call being made is not a valid one")

# A catch all if an AbstractCall's child doesn't define `Base.Expr`
Base.Expr(::AbstractCall) = throw(InvalidCallException)

struct UnaryCall <: AbstractCall
  operator::Union{Symbol, Expr}
  equality::Symbol
  input::Symbol
  output::Symbol
end

# ! WARNING: Do not pass this an inplace function without setting equality to :.=
Base.Expr(c::UnaryCall) = begin
  operator = c.operator
  if c.equality == :.=
    # TODO: Generalize to inplacable functions
    if operator == add_inplace_stub(:⋆₁⁻¹) # Since inverse hodge Geo is a solver
      Expr(:call, c.operator, c.output, c.input)
    elseif operator == :.-
      Expr(c.equality, c.output, Expr(:call, operator, c.input))
    else # TODO: Add check that this operator is a matrix
      Expr(:call, :mul!, c.output, operator, c.input)
    end
  else
    Expr(c.equality, c.output, Expr(:call, operator, c.input))
  end
end

struct BinaryCall <: AbstractCall
  operator::Union{Symbol, Expr}
  equality::Symbol
  input1::Symbol
  input2::Symbol
  output::Symbol
end

# ! WARNING: Do not pass this an inplace function without setting equality to :.=, vice versa
Base.Expr(c::BinaryCall) = begin
  # These operators can be done in-place
  if c.equality == :.= && get_stub(c.operator) == GENSIM_INPLACE_STUB
    return Expr(:call, c.operator, c.output, c.input1, c.input2)
  end
  return Expr(c.equality, c.output, Expr(:call, c.operator, c.input1, c.input2))
end

struct SummationCall <: AbstractCall
  equality::Symbol
  inputs::Vector{Symbol}
  output::Symbol
end

# The output of @code_llvm (.+) of more than 32 variables is inefficient.
Base.Expr(c::SummationCall) = begin
  length(c.inputs) ≤ 32 ?
    Expr(c.equality, c.output, Expr(:call, Expr(:., :+), c.inputs...)) : # (.+)(a,b,c)
    Expr(c.equality, c.output, Expr(:call, :sum, Expr(:vect, c.inputs...))) # sum([a,b,c])
end

struct AllocVecCall <: AbstractCall
  name::Symbol
  form::Symbol
  dimension::Int
  T::DataType
  code_target::AbstractGenerationTarget
end

struct AllocVecCallException <: Exception
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
    _ => throw(AllocVecCallException(c))
  end

  hook_AVC_caching(c, resolved_form, c.code_target)
end

# TODO: Should maybe have default CPU generation be Vector with PreallocTools being opt-in
"""
    hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::CPUBackend)

This hook can be overridden to change the way in which vectors can be preallocated for use by in-place functions.
The AllocVecCall stores the `name` of the vector, the `form` type, the `dimension` of the simulation, the `T` which is the
datatype of the vector, and the `code_target` which is used by multiple dispatch to select a hook.

An example overloaded hook signature would be `hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::UserTarget)`
"""
function hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::CPUBackend)
  :($(Symbol(:__,c.name)) = Decapodes.FixedSizeDiffCache(Vector{$(c.T)}(undef, nparts(mesh, $(QuoteNode(resolved_form))))))
end

# TODO: Allow user to overload these hooks with user-defined code_target
function hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::CUDABackend)
  :($(c.name) = CuVector{$(c.T)}(undef, nparts(mesh, $(QuoteNode(resolved_form)))))
end

"""
    compile_var(alloc_vectors::Vector{AllocVecCall})

This creates the vector allocations that will be used by the simulation body for in-place operations.
"""
function compile_var(alloc_vectors::Vector{AllocVecCall})
  return quote $(Expr.(alloc_vectors)...) end
end

#= function get_form_number(d::SummationDecapode, var_id::Int)
  type = d[var_id, :type]
  if type == :Form0
    return 0
  elseif type == :Form1
    return 1
  elseif type == :Form2
    return 2
  end
  return -1
end

# ! WARNING: This may not work if names are not unique, use above instead
function get_form_number(d::SummationDecapode, var_name::Symbol)
var_id = first(incident(d, var_name, :name))
return get_form_number(d, var_id)
end =#

# TODO: This should be edited when we replace types as symbols with types as Julia types
function is_form(d::SummationDecapode, var_id::Int)
  type = d[var_id, :type]
  return (type == :Form0 || type == :Form1 || type == :Form2 ||
    type == :DualForm0 || type == :DualForm1 || type == :DualForm2)
end

is_form(d::SummationDecapode, var_name::Symbol) = is_form(d, first(incident(d, var_name, :name)))

function getgeneric_type(type::Symbol)
  if (type == :Form0 || type == :Form1 || type == :Form2 ||
    type == :DualForm0 || type == :DualForm1 || type == :DualForm2)
    return :Form
  end
  return type
end

is_literal(d::SummationDecapode, var_id::Int) = (d[var_id, :type] == :Literal)
is_literal(d::SummationDecapode, var_name::Symbol) = is_literal(d, first(incident(d, var_name, :name)))

is_infer(d::SummationDecapode, var_id::Int) = (d[var_id, :type] == :infer)
is_infer(d::SummationDecapode, var_name::Symbol) = is_infer(d, first(incident(d, var_name, :name)))

struct InvalidStubException <: Exception
  stub::Symbol
end

Base.showerror(io::IO, e::InvalidStubException) = print(io, "Stub \"$(e.stub)\" is invalid")

function add_stub(stub_name::Symbol, var_name::Symbol)
  # No empty stubs
  if stub_name == Symbol() || stub_name == NO_STUB_RETURN || !isascii(String(stub_name))
    throw(InvalidStubException(stub_name))
  end

  return Symbol("$(stub_name)_$(var_name)")
end

# ! WARNING: The variable names will be picked up as a stub if no stub was added and the
# ! variable has an "_" itself
function get_stub(var_name::Symbol)
  var_str = String(var_name)
  idx = findfirst("_", var_str)

  if isnothing(idx)
    return NO_STUB_RETURN
  end

  if first(idx) == 1
    throw(InvalidStubException(var_name))
  end

  return Symbol(var_str[begin:first(idx) - 1])
end

add_inplace_stub(var_name::Symbol) = add_stub(GENSIM_INPLACE_STUB, var_name)

const ARITHMETIC_OPS = Set([:+, :*, :-, :/, :.+, :.*, :.-, :./, :^, :.^, :.>, :.<, :.≤, :.≥])

struct InvalidCodeTargetException <: Exception
  code_target::AbstractGenerationTarget
end

Base.showerror(io::IO, e::InvalidCodeTargetException) = print(io, "Provided code target $(e.code_target) is not yet supported in simulations")

"""
    compile_env(d::SummationDecapode, dec_matrices::Vector{Symbol}, con_dec_operators::Set{Symbol}, code_target::AbstractGenerationTarget)

This creates the symbol to function linking for the simulation output. Those run through the `default_dec` backend
expect both an in-place and an out-of-place variant in that order. User defined operations only support out-of-place.
"""
function compile_env(d::SummationDecapode, dec_matrices::Vector{Symbol}, con_dec_operators::Set{Symbol}, code_target::AbstractGenerationTarget)
  defined_ops = deepcopy(con_dec_operators)

  defs = quote end

  for op in dec_matrices
    op in defined_ops && continue

    quote_op = QuoteNode(op)
    mat_op = add_stub(GENSIM_INPLACE_STUB, op)

    # TODO: Add support for user-defined code targets
    default_generation = @match code_target begin
      ::CPUBackend => :default_dec_matrix_generate
      ::CUDABackend => :default_dec_cu_matrix_generate
      _ => throw(InvalidCodeTargetException(code_target))
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

struct AmbiguousNameException <: Exception
  name::Symbol
  indices::Vector{Int}
end

Base.showerror(io::IO, e::AmbiguousNameException) = begin
  if isempty(e.indices)
    print(io, "Name \"$(e.name)\" is does not exist")
  else
    print(io, "Name \"$(e.name)\" is repeated at indices $(e.indices) and is ambiguous")
  end
end

struct InvalidDecaTypeException <: Exception
  name::Symbol
  type::Symbol
end

Base.showerror(io::IO, e::InvalidDecaTypeException) = print(io, "Variable \"$(e.name)\" has invalid type \"$(e.type)\"")

"""
    get_vars_code(d::SummationDecapode, vars::Vector{Symbol}, ::Type{stateeltype}, code_target::AbstractGenerationTarget) where stateeltype

This initalizes all input variables according to their Decapodes type.
"""
function get_vars_code(d::SummationDecapode, vars::Vector{Symbol}, ::Type{stateeltype}, code_target::AbstractGenerationTarget) where stateeltype
  stmts = quote end

  foreach(vars) do s
    # If name is not unique (or not just literals) then error
    found_names_idxs = incident(d, s, :name)

    # TODO: we should handle the case of same literals better
    is_singular = length(found_names_idxs) == 1
    is_all_literals = all(d[found_names_idxs, :type] .== :Literal)

    if isempty(found_names_idxs) || (!is_singular && !is_all_literals)
      throw(AmbiguousNameException(s, found_names_idxs))
    end

    s_type = is_all_literals ? :Literal : getgeneric_type(d[only(found_names_idxs), :type])

    # Literals don't need assignments, because they are literals, but we stored them as Symbols.
    # TODO: we should fix that upstream so that we don't need this.
    line = @match s_type begin
      :Literal => :($s = $(parse(stateeltype, String(s))))
      :Constant => :($s = p.$s)
      :Parameter => :($s = (p.$s)(t))
      _ => hook_GVC_get_form(s, s_type, code_target) # ! WARNING: This assumes a form
      # _ => throw(InvalidDecaTypeException(s, s_type)) # TODO: Use this for invalid types
    end
    push!(stmts.args, line)
  end
  return stmts
end

# TODO: Expand on this to be able to handle vector and ComponentArrays inputs
function hook_GVC_get_form(var_name::Symbol, var_type::Symbol, ::Union{CPUBackend, CUDABackend})
  return :($var_name = u.$var_name)
end

"""
    set_tanvars_code(d::SummationDecapode)

This function creates the code that sets the value of the Tvars at the end of the code
"""
function set_tanvars_code(d::SummationDecapode, code_target::AbstractGenerationTarget)
  stmts = quote end

  tanvars = [(d[e, [:src,:name]], d[e, [:tgt,:name]]) for e in incident(d, :∂ₜ, :op1)]
  foreach(tanvars) do (s,t)
    push!(stmts.args, hook_STC_settvar(s, t, code_target))
  end
  return stmts
end

# TODO: Why do we need a QuoteNode for this?
"""
    hook_STC_settvar(src_name::Symbol, tgt_name::Symbol, ::Union{CPUBackend, CUDABackend})

This hook is meant to control how data is set into the tangent variables after a simulation function
execution. It expects the `state_name`, the name of the original state variable, the `tgt_name` which
is the name of the variable whose data will be stored and a code target.
"""
function hook_STC_settvar(state_name::Symbol, tgt_name::Symbol, ::Union{CPUBackend, CUDABackend})
  ssymb = QuoteNode(state_name)
  return :(setproperty!(du, $ssymb, $tgt_name))
end

const PROMOTE_ARITHMETIC_MAP = Dict(:(+) => :.+,
                                    :(-) => :.-,
                                    :(*) => :.*,
                                    :(/) => :./,
                                    :(^) => :.^,
                                    :(=) => :.=,
                                    :.+ => :.+,
                                    :.- => :.-,
                                    :.* => :.*,
                                    :./ => :./,
                                    :.^ => :.^,
                                    :.= => :.=)

"""
    compile(d::SummationDecapode, inputs::Vector{Symbol}, alloc_vectors::Vector{AllocVecCall}, optimizable_dec_operators::Set{Symbol}, dimension::Int, stateeltype::DataType, code_target::AbstractGenerationTarget, preallocate::Bool)

Function that compiles the computation body. `d` is the input Decapode, `inputs` is a vector of state variables and literals,
`alloc_vec` should be empty when passed in, `optimizable_dec_operators` is a collection of all DEC operator symbols that can use special
in-place methods, `dimension` is the dimension of the problem (usually 1 or 2), `stateeltype` is the type of the state elements
(usually Float32 or Float64), `code_target` determines what architecture the code is compiled for (either CPU or CUDA), and `preallocate`
which is set to `true` by default and determines if intermediate results can be preallocated..
"""
function compile(d::SummationDecapode, inputs::Vector{Symbol}, alloc_vectors::Vector{AllocVecCall}, optimizable_dec_operators::Set{Symbol}, dimension::Int, stateeltype::DataType, code_target::AbstractGenerationTarget, preallocate::Bool)
  # Get the Vars of the inputs (probably state Vars).
  visited_Var = falses(nparts(d, :Var))

  # TODO: Pass in state indices instead of names
  input_numbers = reduce(vcat, incident(d, inputs, :name))

  visited_Var[input_numbers] .= true
  visited_Var[incident(d, :Literal, :type)] .= true

  # TODO: Collect these visited arrays into one structure indexed by :Op1, :Op2, and :Σ
  visited_1 = falses(nparts(d, :Op1))
  visited_2 = falses(nparts(d, :Op2))
  visited_Σ = falses(nparts(d, :Σ))

  # FIXME: this is a quadratic implementation of topological_sort inlined in here.
  op_order = AbstractCall[]

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
        if preallocate && is_form(d, t)
          if operator in optimizable_dec_operators
            equality = PROMOTE_ARITHMETIC_MAP[equality]
            operator = add_stub(GENSIM_INPLACE_STUB, operator)
            push!(alloc_vectors, AllocVecCall(tname, d[t, :type], dimension, stateeltype, code_target))

          elseif operator == :(-) || operator == :.-
            equality = PROMOTE_ARITHMETIC_MAP[equality]
            operator = PROMOTE_ARITHMETIC_MAP[operator]
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

        # TODO: Check to make sure that this logic never breaks
        if preallocate && is_form(d, r)
          if operator == :(+) || operator == :(-) || operator == :.+ || operator == :.-
            operator = PROMOTE_ARITHMETIC_MAP[operator]
            equality = PROMOTE_ARITHMETIC_MAP[equality]
            push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))

          # TODO: Do we want to support the ability of a user to use the backslash operator?
          elseif operator == :(*) || operator == :(/) || operator == :.* || operator == :./
            # ! WARNING: This part may break if we add more compiler types that have different
            # ! operations for basic and broadcast modes, e.g. matrix multiplication vs broadcast
            if !is_infer(d, arg1) && !is_infer(d, arg2)
              operator = PROMOTE_ARITHMETIC_MAP[operator]
              equality = PROMOTE_ARITHMETIC_MAP[equality]
              push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))
            end
          elseif operator in optimizable_dec_operators
            operator = add_stub(GENSIM_INPLACE_STUB, operator)
            equality = PROMOTE_ARITHMETIC_MAP[equality]
            push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))
          end
        end

        # TODO: Clean this in another PR (with a @match maybe).
        if operator == :(*)
          operator = PROMOTE_ARITHMETIC_MAP[operator]
        end
        if operator == :(-)
          operator = PROMOTE_ARITHMETIC_MAP[operator]
        end
        if operator == :(/)
          operator = PROMOTE_ARITHMETIC_MAP[operator]
        end
        if operator == :(^)
          operator = PROMOTE_ARITHMETIC_MAP[operator]
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

        # operator = :(+)
        operator = :.+
        equality = :(=)

        # If result is a known form, broadcast addition
        if preallocate && is_form(d, r)
          operator = PROMOTE_ARITHMETIC_MAP[operator]
          equality = PROMOTE_ARITHMETIC_MAP[equality]
          push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype, code_target))
        end

        visited_Σ[op] = true
        visited_Var[r] = true
        c = SummationCall(equality, argnames, rname)
        push!(op_order, c)
      end
    end
  end

  eq_exprs = Expr.(op_order)
end

"""
    post_process_vector_allocs(alloc_vecs::Vector{AllocVecCall}, code_target::AbstractGenerationTarget)

This deals with any post processing needed by the allocations, like if data needs to be retrieved
from a special cache.
"""
function post_process_vector_allocs(alloc_vecs::Vector{AllocVecCall}, code_target::AbstractGenerationTarget)
  list_exprs = Expr[]
  foreach(alloc_vecs) do alloc_vec
    hook_PPVA_data_handle!(list_exprs, alloc_vec, code_target)
  end

  cache_exprs = quote end
  append!(cache_exprs.args, list_exprs)
  return cache_exprs
end

"""
    hook_PPVA_data_handle!(cache_exprs::Vector{Expr}, alloc_vec::AllocVecCall, ::CPUBackend)

This hook determines if preallocated vectors need to be be handled in a special manner everytime
before a function run. This is useful in the example of using `FixedSizeDiffCache` from `PreallocationTools.jl`.

This hook is passed in `cache_exprs` which is the collection of exprs to be pasted, `alloc_vec` which is an
`AllocVecCall` that stores information about the allocated vector and a code target.
"""
function hook_PPVA_data_handle!(cache_exprs::Vector{Expr}, alloc_vec::AllocVecCall, ::CPUBackend)
  line = :($(alloc_vec.name) = (Decapodes.get_tmp($(Symbol(:__,alloc_vec.name)), u)))
  push!(cache_exprs, line)
end

function hook_PPVA_data_handle!(cache_exprs::Vector{Expr}, alloc_vec::AllocVecCall, ::CUDABackend)
  return
end

"""
    resolve_types_compiler!(d::SummationDecapode)

Converts `Constant` and `Parameter` types to `infer` since this is essentially what they are
to the compiler.
"""
function resolve_types_compiler!(d::SummationDecapode)
  d[:type] = map(d[:type]) do x
    if x == :Constant || x == :Parameter
      return :infer
    end
    return x
  end
end

"""
    replace_names_compiler!(d::SummationDecapode)

This makes easy function name conversions in the Decapode
"""
function replace_names_compiler!(d::SummationDecapode)
  dec_op1 = Pair{Symbol, Any}[]
  dec_op2 = Pair{Symbol, Symbol}[(:∧₀₀ => :.*)]
  replace_names!(d, dec_op1, dec_op2)
end

# TODO: This should be extended to accept user rules
"""
    infer_overload_compiler!(d::SummationDecapode, dimension::Int)

A combined `infer_types` and `resolve_overloads` pipeline with default DEC rules.
"""
function infer_overload_compiler!(d::SummationDecapode, dimension::Int)
  if dimension == 1
    infer_types!(d, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(d, op1_res_rules_1D, op2_res_rules_1D)
  elseif dimension == 2
    infer_types!(d, op1_inf_rules_2D, op2_inf_rules_2D)
    resolve_overloads!(d, op1_res_rules_2D, op2_res_rules_2D)
  end
end

"""
    init_dec_matrices!(d::SummationDecapode, dec_matrices::Vector{Symbol}, optimizable_dec_operators::Set{Symbol})

Collects all DEC operators that are concrete matrices.
"""
function init_dec_matrices!(d::SummationDecapode, dec_matrices::Vector{Symbol}, optimizable_dec_operators::Set{Symbol})
  for op_name in vcat(d[:op1], d[:op2])
    if op_name in optimizable_dec_operators
      push!(dec_matrices, op_name)
    end
  end
end

"""
    link_contract_operators(d::SummationDecapode, con_dec_operators::Set{Symbol}, stateeltype::DataType, code_target::AbstractGenerationTarget)

Collects arrays of DEC matrices together, replaces the array with a generated function name and computes the contracted multiplication
"""
function link_contract_operators(d::SummationDecapode, con_dec_operators::Set{Symbol}, stateeltype::DataType, code_target::AbstractGenerationTarget)

  contract_defs = quote end

  compute_to_name = Dict()
  curr_id = 1

  for op1_id in parts(d, :Op1)
    op1_name = d[op1_id, :op1]
    if isa(op1_name, AbstractArray)
      computation = reverse!(map(x -> add_inplace_stub(x), op1_name))
      compute_key = join(computation, " * ")

      computation_name = get(compute_to_name, compute_key, :Error)
      if computation_name == :Error
        computation_name = add_stub(Symbol("GenSim-ConMat"), Symbol(curr_id))
        get!(compute_to_name, compute_key, computation_name)
        push!(con_dec_operators, computation_name)

        expr_line = hook_LCO_inplace(computation_name, computation, stateeltype, code_target)
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
function hook_LCO_inplace(computation_name::Symbol, computation::Vector{Symbol}, float_type::DataType, ::CPUBackend)
  return :($(add_inplace_stub(computation_name)) = $(Expr(:call, :*, computation...)))
end

function generate_parentheses_multiply(list)
  if length(list) == 1
      return list[1]
  else
      return Expr(:call, :*, generate_parentheses_multiply(list[1:end-1]), list[end])
  end
end

function hook_LCO_inplace(computation_name::Symbol, computation::Vector{Symbol}, float_type::DataType, ::CUDABackend)
  return :($(add_inplace_stub(computation_name)) = $(generate_parentheses_multiply(computation)))
end

struct UnsupportedDimensionException <: Exception
  dim::Int
end

Base.showerror(io::IO, e::UnsupportedDimensionException) = print(io, "Decapodes does not support dimension $(e.dim) simulations")

struct UnsupportedStateeltypeException <: Exception
  type::DataType
end

Base.showerror(io::IO, e::UnsupportedStateeltypeException) = print(io, "Decapodes does not support state element types as $(e.type), only Float32 or Float64")

"""
    gensim(user_d::SummationDecapode, input_vars::Vector{Symbol}; dimension::Int=2, stateeltype::DataType = Float64, code_target::AbstractGenerationTarget = CPUTarget(), preallocate::Bool = true)

Generates the entire code body for the simulation function. The returned simulation function can then be combined with a mesh, provided by `CombinatorialSpaces`, and a function describing symbol 
to operator mappings to return a simulator that can be used to solve the represented equations given initial conditions.
  
**Arguments:**
  
`user_d`: The user passed Decapode for which simulation code will be generated. (This is not modified) 

`input_vars` is the collection of variables whose values are known at the beginning of the simulation. (Defaults to all state variables and literals in the Decapode)

**Keyword arguments:**

`dimension`: The dimension of the problem. (Defaults to `2`)(Must be `1` or `2`)

`stateeltype`: The element type of the state forms. (Defaults to `Float64`)(Must be `Float32` or `Float64`)

`code_target`: The intended architecture target for the generated code. (Defaults to `CPUTarget()`)(Use `CUDATarget()` for NVIDIA CUDA GPUs)

`preallocate`: Enables(`true`)/disables(`false`) pre-allocated caches for intermediate computations. Some functions, such as those that determine Jacobian sparsity patterns, or perform auto-differentiation, may require this to be disabled. (Defaults to `true`)
"""
function gensim(user_d::SummationDecapode, input_vars::Vector{Symbol}; dimension::Int=2, stateeltype::DataType = Float64, code_target::AbstractGenerationTarget = CPUTarget(), preallocate::Bool = true)

  (dimension == 1 || dimension == 2) ||
    throw(UnsupportedDimensionException(dimension))

  (stateeltype == Float32 || stateeltype == Float64) ||
    throw(UnsupportedStateeltypeException(stateeltype))

  # Explicit copy for safety
  gen_d = deepcopy(user_d)

  recognize_types(gen_d)

  # Makes copy
  gen_d = expand_operators(gen_d)

  dec_matrices = Vector{Symbol}()
  alloc_vectors = Vector{AllocVecCall}()

  vars = get_vars_code(gen_d, input_vars, stateeltype, code_target)
  tars = set_tanvars_code(gen_d, code_target)

  # We need to run this after we grab the constants and parameters out
  infer_overload_compiler!(gen_d, dimension)
  resolve_types_compiler!(gen_d)
  infer_overload_compiler!(gen_d, dimension)

  # This should probably be followed by an expand_operators
  replace_names_compiler!(gen_d)
  open_operators!(gen_d, dimension = dimension)
  infer_overload_compiler!(gen_d, dimension)

  # This will generate all of the fundemental DEC operators present
  optimizable_dec_operators = Set([:⋆₀, :⋆₁, :⋆₂, :⋆₀⁻¹, :⋆₂⁻¹,
                                  :d₀, :d₁, :dual_d₀, :d̃₀, :dual_d₁, :d̃₁,
                                  :avg₀₁])
  extra_dec_operators = Set([:⋆₁⁻¹, :∧₀₁, :∧₁₀, :∧₁₁, :∧₀₂, :∧₂₀])

  init_dec_matrices!(gen_d, dec_matrices, union(optimizable_dec_operators, extra_dec_operators))

  # This contracts matrices together into a single matrix
  contracted_dec_operators = Set{Symbol}()
  contract_operators!(gen_d, allowable_ops = optimizable_dec_operators)
  cont_defs = link_contract_operators(gen_d, contracted_dec_operators, stateeltype, code_target)

  union!(optimizable_dec_operators, contracted_dec_operators, extra_dec_operators)

  # Compilation of the simulation
  equations = compile(gen_d, input_vars, alloc_vectors, optimizable_dec_operators, dimension, stateeltype, code_target, preallocate)
  data = post_process_vector_allocs(alloc_vectors, code_target)

  func_defs = compile_env(gen_d, dec_matrices, contracted_dec_operators, code_target)
  vect_defs = compile_var(alloc_vectors)

  quote
    (mesh, operators, hodge=GeometricHodge()) -> begin
      $func_defs
      $cont_defs
      $vect_defs
      f(du, u, p, t) = begin
        $vars
        $data
        $(equations...)
        $tars
        return nothing
      end;
    end
  end
end

gather_inputs(d::SummationDecapode) = vcat(infer_state_names(d), d[incident(d, :Literal, :type), :name])

gensim(c::Collage; dimension::Int=2) =
gensim(collate(c); dimension=dimension)

gensim(d::SummationDecapode; kwargs...) =
  gensim(d, gather_inputs(d); kwargs...)

evalsim(args...; kwargs...) = eval(gensim(args...; kwargs...))

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
