using CombinatorialSpaces
using ComponentArrays
using ForwardDiff
using LinearAlgebra
using MLStyle
using PreallocationTools

const GENSIM_INPLACE_STUB = Symbol("GenSim-M")
const NO_STUB_RETURN = Symbol("NOSTUB")

abstract type AbstractGenerationTarget end

abstract type CPUBackend <: AbstractGenerationTarget end
abstract type CUDABackend <: AbstractGenerationTarget end
# TODO: Test that AbstractGenerationTargets are user-extendable.

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

opt_generator_function(code_target::AbstractGenerationTarget) = throw(InvalidCodeTargetException(code_target))
opt_generator_function(::CPUBackend) = :default_dec_matrix_generate
opt_generator_function(::CUDABackend) = :default_dec_cu_matrix_generate

generator_function(code_target::AbstractGenerationTarget) = throw(InvalidCodeTargetException(code_target))
generator_function(::CPUBackend) = :default_dec_generate
generator_function(::CUDABackend) = :default_dec_cu_generate

# TODO: This function should be handled with dispatch.
"""    compile_env(d::SummationDecapode, present_dec_ops::Vector{Symbol}, contracted_ops::Vector{Symbol}, code_target::AbstractGenerationTarget)

Emit code to define functions given operator Symbols.

Default operations return a tuple of an in-place and an out-of-place function. User-defined operations return an out-of-place function.
"""
function compile_env(d::SummationDecapode, present_dec_ops::Set{Symbol}, contracted_ops::Vector{Symbol}, code_target::AbstractGenerationTarget)

  defs = quote end

  all_ops = d[:op1] ∪ d[:op2] ∪ present_dec_ops
  avoid_ops = contracted_ops ∪ [DerivOp] ∪ ARITHMETIC_OPS

  for op in setdiff(all_ops, avoid_ops)
    quote_op = QuoteNode(op)
    def = @match op begin
      if op in optimizable(code_target) end => :(($(add_inplace_stub(op)), $op) = $(opt_generator_function(code_target))(mesh, $quote_op, hodge))
      if op in non_optimizable(code_target) end => :($op = $(generator_function(code_target))(mesh, $quote_op, hodge))
      _ => :($op = operators(mesh, $quote_op))
    end
    push!(defs.args, def)
  end

  defs
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
      :Constant => :($s = __p__.$s)
      :Parameter => :($s = (__p__.$s)(__t__))
      _ => hook_GVC_get_form(s, s_type, code_target) # ! WARNING: This assumes a form
      # _ => throw(InvalidDecaTypeException(s, s_type)) # TODO: Use this for invalid types
    end
    push!(stmts.args, line)
  end
  return stmts
end

# TODO: Expand on this to be able to handle vector and ComponentArrays inputs
function hook_GVC_get_form(var_name::Symbol, var_type::Symbol, ::Union{CPUBackend, CUDABackend})
  return :($var_name = __u__.$var_name)
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
  return :(setproperty!(__du__, $ssymb, $tgt_name))
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
    compile(d::SummationDecapode, inputs::Vector{Symbol}, inplace_dec_ops::Set{Symbol}, dimension::Int, stateeltype::DataType, code_target::AbstractGenerationTarget, preallocate::Bool)

Function that compiles the computation body. `d` is the input Decapode, `inputs` is a vector of state variables and literals,
`inplace_dec_ops` is a collection of all DEC operator symbols that can use special
in-place methods, `dimension` is the dimension of the problem (usually 1 or 2), `stateeltype` is the type of the state elements
(usually Float32 or Float64), `code_target` determines what architecture the code is compiled for (either CPU or CUDA), and `preallocate`
which is set to `true` by default and determines if intermediate results can be preallocated..
"""
function compile(d::SummationDecapode, inputs::Vector{Symbol}, inplace_dec_ops::Set{Symbol}, dimension::Int, stateeltype::DataType, code_target::AbstractGenerationTarget, preallocate::Bool)
  alloc_vectors = Vector{AllocVecCall}()
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

        if preallocate && is_form(d, t)
          if operator in inplace_dec_ops
            equality = PROMOTE_ARITHMETIC_MAP[equality]
            operator = add_inplace_stub(operator)
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
          elseif operator in inplace_dec_ops
            operator = add_inplace_stub(operator)
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

  Expr.(op_order), alloc_vectors
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
  line = :($(alloc_vec.name) = (Decapodes.get_tmp($(Symbol(:__,alloc_vec.name)), __u__)))
  push!(cache_exprs, line)
end

function hook_PPVA_data_handle!(cache_exprs::Vector{Expr}, alloc_vec::AllocVecCall, ::CUDABackend)
  return
end

"""
    convert_cs_ps_to_infer!(d::SummationDecapode)

Convert `Constant` and `Parameter` types to `infer`.
"""
function convert_cs_ps_to_infer!(d::SummationDecapode)
  cs_ps = incident(d, [:Constant, :Parameter], :type)
  d[vcat(cs_ps...), :type] = :infer
end

# TODO: This should be extended to accept user rules
"""
    infer_overload_compiler!(d::SummationDecapode, dimension::Int)

A combined `infer_types` and `resolve_overloads` pipeline with default DEC rules.
"""
function infer_overload_compiler!(d::SummationDecapode, dimension::Int)
  infer_types!(d, dim = dimension)
  resolve_overloads!(d, dim = dimension)
end

"""    link_contracted_operators!(d::SummationDecapode, code_target::AbstractGenerationTarget)

Emit code to pre-multiply unique sequences of matrix operations, and rename corresponding operations.
"""
function link_contracted_operators!(d::SummationDecapode, code_target::AbstractGenerationTarget)
  contracted_defs = quote end
  contracted_ops = Symbol[]
  chain_idxs = findall(x -> x isa AbstractArray, d[:op1])

  for (i, chain) in enumerate(unique(d[chain_idxs, :op1]))
    LHS = add_stub(Symbol("GenSim-ConMat"), Symbol(i-1))
    RHS = reverse!(add_inplace_stub.(chain))

    push!(contracted_ops, LHS)
    push!(contracted_defs.args, mat_def_expr(LHS, RHS, code_target), mat_mul_func_expr(LHS))
    d[findall(==(chain), d[:op1]), :op1] = LHS
  end
  contracted_defs, contracted_ops
end

# Given the name of a matrix, return an Expr that multiplies by that matrix.
mat_mul_func_expr(mat_name) =
  :($mat_name = x -> $(add_inplace_stub(mat_name)) * x)

# Given the name and factors of a matrix, return an Expr that defines that matrix.
mat_def_expr(computation_name::Symbol, factors::Vector{Symbol}, ::CPUBackend) =
  :($(add_inplace_stub(computation_name)) = *($(factors...)))

nested_mul(factors) =
  length(factors) == 1 ?
    factors[begin] :
    Expr(:call, :*, nested_mul(factors[begin:end-1]), factors[end])

mat_def_expr(computation_name::Symbol, factors::Vector{Symbol}, ::CUDABackend) =
  :($(add_inplace_stub(computation_name)) = $(nested_mul(factors)))

struct UnsupportedDimensionException <: Exception
  dim::Int
end

Base.showerror(io::IO, e::UnsupportedDimensionException) = print(io, "Decapodes does not support dimension $(e.dim) simulations")

struct UnsupportedStateeltypeException <: Exception
  type::DataType
end

Base.showerror(io::IO, e::UnsupportedStateeltypeException) = print(io, "Decapodes does not support state element types as $(e.type), only Float32 or Float64")

const MATRIX_OPTIMIZABLE_DEC_OPERATORS = Set([:⋆₀, :⋆₁, :⋆₂, :⋆₀⁻¹, :⋆₂⁻¹,
                                              :d₀, :d₁, :dual_d₀, :d̃₀, :dual_d₁, :d̃₁,
                                              :avg₀₁, :♭♯])

const NONMATRIX_OPTIMIZABLE_DEC_OPERATORS = Set([:⋆₁⁻¹, :∧₀₁, :∧₁₀, :∧₁₁, :∧₀₂, :∧₂₀])


const NON_OPTIMIZABLE_CPU_OPERATORS = Set([:♯ᵖᵈ, :♯ᵖᵖ, :♯ᵈᵈ, :♭ᵈᵖ,
                                           :∧ᵖᵈ₁₁, :∧ᵖᵈ₀₁, :∧ᵈᵖ₁₁, :∧ᵈᵖ₁₀, :∧ᵈᵈ₁₁, :∧ᵈᵈ₁₀, :∧ᵈᵈ₀₁,
                                           :ι₁₁, :ι₁₂, :ℒ₁, :Δᵈ₀ , :Δᵈ₁, :Δ₀⁻¹, :neg, :mag, :norm])
const NON_OPTIMIZABLE_CUDA_OPERATORS = Set{Symbol}()

const DEC_GEN_OPTIMIZABLE_OPERATORS = MATRIX_OPTIMIZABLE_DEC_OPERATORS ∪ NONMATRIX_OPTIMIZABLE_DEC_OPERATORS

optimizable(code_target::AbstractGenerationTarget) = throw(InvalidCodeTargetException(code_target))
optimizable(::CPUBackend) = DEC_GEN_OPTIMIZABLE_OPERATORS
optimizable(::CUDABackend) = DEC_GEN_OPTIMIZABLE_OPERATORS

non_optimizable(code_target::AbstractGenerationTarget) = throw(InvalidCodeTargetException(code_target))
non_optimizable(::CPUBackend) = NON_OPTIMIZABLE_CPU_OPERATORS
non_optimizable(::CUDABackend) = NON_OPTIMIZABLE_CUDA_OPERATORS

dec_operator_set(code_target::AbstractGenerationTarget) = optimizable(code_target) ∪ non_optimizable(code_target)

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

`contract`: Enables(`true`)/disables(`false`) pre-computation of matrix-matrix multiplications for chains of such operators. This feature can interfere with certain auto-differentiation methods, in which case this can be disabled. (Defaults to `true`)

`multigrid`: Enables multigrid methods during code generation. If `true`, then the function produced by `gensim` will expect a `PrimalGeometricMapSeries`. (Defaults to `false`)
"""
function gensim(user_d::SummationDecapode, input_vars::Vector{Symbol}; dimension::Int=2, stateeltype::DataType = Float64, code_target::AbstractGenerationTarget = CPUTarget(), preallocate::Bool = true, contract::Bool = true, multigrid::Bool = false)
  (dimension == 1 || dimension == 2) ||
    throw(UnsupportedDimensionException(dimension))

  (stateeltype == Float32 || stateeltype == Float64) ||
    throw(UnsupportedStateeltypeException(stateeltype))

  # Explicit copy for safety
  d = deepcopy(user_d)

  recognize_types(d)

  # Makes copy
  d = expand_operators(d)

  vars = get_vars_code(d, input_vars, stateeltype, code_target)
  tars = set_tanvars_code(d, code_target)

  infer_overload_compiler!(d, dimension)
  convert_cs_ps_to_infer!(d)
  infer_overload_compiler!(d, dimension)

  # XXX: expand_operators should be called if any replacement is a chain of operations.
  replace_names!(d, Pair{Symbol, Any}[], Pair{Symbol, Symbol}[(:∧₀₀ => :.*)])
  open_operators!(d, dimension = dimension)
  infer_overload_compiler!(d, dimension)

  present_dec_ops = Set{Symbol}(dec_operator_set(code_target) ∩ (d[:op1] ∪ d[:op2]))

  # This contracts matrices together into a single matrix
  contract && contract_operators!(d, white_list = MATRIX_OPTIMIZABLE_DEC_OPERATORS)
  contracted_defs, contracted_ops = link_contracted_operators!(d, code_target)

  # Combination of already in-place dec operators and newly contracted matrices
  inplace_dec_ops = union(optimizable(code_target), contracted_ops)

  # Compilation of the simulation
  equations, alloc_vectors = compile(d, input_vars, inplace_dec_ops, dimension, stateeltype, code_target, preallocate)
  data = post_process_vector_allocs(alloc_vectors, code_target)

  func_defs = compile_env(d, present_dec_ops, contracted_ops, code_target)
  vect_defs = quote $(Expr.(alloc_vectors)...) end

  multigrid_defs = quote end
  multigrid && push!(multigrid_defs.args, :(mesh = finest_mesh(mesh)))

  quote
    (mesh, operators, hodge=GeometricHodge()) -> begin
      $func_defs
      $contracted_defs
      $multigrid_defs
      $vect_defs
      f(__du__, __u__, __p__, __t__) = begin
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
