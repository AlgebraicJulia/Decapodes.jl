using CombinatorialSpaces
using ComponentArrays
using OrdinaryDiffEq
using GeometryBasics
using LinearAlgebra
using Base.Iterators
using Catlab
using MLStyle
import Catlab.Programs.GenerateJuliaPrograms: compile

using Parameters # for default values in structs

@enum Stub in_place_stub=Symbol("GenSim-M") contract_stub=Symbol("GenSim-ConMat")
@enum Operation addition=:A multiplication=:M

abstract type AbstractOperator end

mutable struct UnaryOperator <: AbstractOperator
  op::Symbol
  stubs::Vector{AbstractStub}
  is_optimizable::Bool
end

#= JAMES
@data OpType begin
  Addition
  Mult
end

@data Operator <: AbstractOperator begin
  Unary(op::Symbol, broadcasted::Bool, optimizable::Bool)
  Binary(op::Symbol, optype::OpType, broadcasted::Bool, optimizable::Bool)
end
=#

mutable struct BinaryOperator <: AbstractOperator
  op::Symbol
  stubs::Vector{AbstractStub}
  optype::Operation # addition or multiplication (and their inverses)
  is_negation::Bool
  is_broadcasted::Bool
  is_optimizable::Bool
end

abstract type AbstractCall end

# TODO: Do not use abstract types as params.
struct UnaryCall <: AbstractCall
  operator::AbstractOperator
  equality::Any
  input::Any
  output::Any
end

function to_op(d::SummationDecapode, op::Int, ::Type{Val{:op1}})
  UnaryOperator(d[op, :op1], [], false)
end

function to_op(d::SummationDecapode, op::Int, ::Type{Val{:op2}})
  @match d[op, :op2] begin
    op::Symbol && if op ∈ [:+, :-] => BinaryOperator(op, [], :A, false, false)
    op::Symbol && if op ∈ [:.+, :.-] => BinaryOperator(op, [], :A, false, true)
    op::Symbol && if op ∈ [:*, :/] => BinaryOperator(op, [], :M, false, false)
    op::Symbol && if op ∈ [:.*, :./] => BinaryOperator(op, [], :M, false, true)
    _ => nothing
  end
end


# per MLStyle manual, we need to use as_record in order to match
@as_record UnaryCall
@as_record UnaryOperator

# TODO: Add back support for contract operators
Base.Expr(call::UnaryCall) = begin
  @match call begin
    UnaryCall(; operator=UnOp(:⋆₁⁻¹, [in_place_stub]), equality=:.=) => Expr(:call, call.operator, c.output, c.input)
    UnaryCall(; equality=:.=) => Expr(:call, :mul!, c.output, c.operator, c.input)
    _ => Expr(call.equality, call.output, Expr(:call, c.operator, c.input))
  end
end

struct BinaryCall <: AbstractCall
  operator::AbstractOperator
  equality::Any
  input1::Any
  input2::Any
  output::Any
end

@as_record BinaryCall

# TODO: After getting rid of AppCirc2, do we need this check?
Base.Expr(call::BinaryCall) = begin
  @match call begin
    BinaryCall(; operator=BinaryOperator(; stubs=[in_place_stub]), equality=:.=) => 
      Expr(:call, call.operator, c.output, call.input1, c.input2)
    _ => Expr(call.equality, call.output, Expr(:call, call.operator, call.input1, call.input2))
  end
end

struct VarargsCall <: AbstractCall
  operator::Any
  equality::Any
  inputs::Any
  output::Any
end

Base.Expr(c::VarargsCall) = begin
  return Expr(c.equality, c.output, Expr(:call, c.operator, c.inputs...))
end

struct AllocVecCall <: AbstractCall
  name::Any
  form::Any
  dimension::Any
  T::Any
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

  :(
    $(Symbol(:__, c.name)) = Decapodes.FixedSizeDiffCache(
      Vector{$(c.T)}(undef, nparts(mesh, $(QuoteNode(resolved_form)))),
    )
  )
end

# WARNING: This may not work if names are not unique, use above instead
function get_form_number(d::SummationDecapode, var_name::Symbol)
var_id = first(incident(d, var_name, :name))
return get_form_number(d, var_id)
end =#

struct Form
  dim::Int
  is_dual::Bool
end

forms = [Form(x,y) for x ∈ 1:3 for y ∈ 0:1]

function is_form(d::SummationDecapode, var_id::Int)
  return d[var_id, :type] ∈ [:Form0, :Form1, :Form2, :DualForm0, :DualForm1, :DualForm2]
end

is_form(d::SummationDecapode, var_name::Symbol) =
  is_form(d, first(incident(d, var_name, :name)))

is_literal(d::SummationDecapode, var_id::Int) = (d[var_id, :type] == :Literal)
is_literal(d::SummationDecapode, var_name::Symbol) =
  is_literal(d, first(incident(d, var_name, :name)))

is_infer(d::SummationDecapode, var_id::Int) = (d[var_id, :type] == :infer)
is_infer(d::SummationDecapode, var_name::Symbol) =
  is_infer(d, first(incident(d, var_name, :name)))

# TODO we should just be appending symbols to the stubs field 
add_stub(stub_name::Symbol, op_name::Symbol) = return Symbol("$(stub_name)_$(op_name)")
add_stub(op_name::Symbol, Type{Val{in_place_stub}}) = add_stub(in_place_stub, op_name)
add_stub(op_name::Symbol, Type{Val{contract_stub}}) = add_stub(contract_stub, op_name)

function get_stub(var_name::Symbol)
  var_str = String(var_name)
  idx = findfirst("_", var_str)
  if (isnothing(idx) || first(idx) == 1)
    return nothing
  end
  return Symbol(var_str[begin:first(idx)-1])
end

# This will be the function and matrix generation
function compile_env(
  d::AbstractNamedDecapode,
  dec_matrices::Vector{Symbol},
  con_dec_operators::Set{Symbol},
)
  # TODO 
  assumed_ops = Set([:+, :*, :-, :/, :.+, :.*, :.-, :./])
  defined_ops = Set()

  defs = quote end

  # TODO can we just use symdiff?
  for op in dec_matrices
    if (op in defined_ops)
      continue
    end

    quote_op = QuoteNode(op)
    mat_op = add_stub(in_place_stub, op)
    # Could turn this into a special call
    def = :(($mat_op, $op) = default_dec_matrix_generate(mesh, $quote_op, hodge))
    push!(defs.args, def)

    push!(defined_ops, op)
  end
  
  # LUKE
  op1s = map(filter(x -> x \notin [con_dec_operators…, DerivOp], unique(d[:op1]))) do op
    :($op = operators(mesh, $(QuoteNode(op))))
  end
  append!(defs.args, op1s)

  # TODO can we just use setdiff?
  for op in setdiff(d[:op1], con_dec_operators, defined_ops, [DerivOp])
    # if op == DerivOp
    #   continue
    # end
    #
    # if (op in con_dec_operators || op in defined_ops)
    #   continue
    # end

    ops = QuoteNode(op)
    def = :($op = operators(mesh, $ops))

    push!(defs.args, def)

    push!(defined_ops, op)
  end
  for op in setdiff(d[:op2], assumed_ops, defined_ops)
    # TODO setdiff?
    if op in assumed_ops || op in defined_ops
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
  return quote
    $(Expr.(alloc_vectors)...)
  end
end

# This is the block of parameter setting inside f
# TODO: Pass this an extra type parameter that sets the size of the Floats
get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol}) =
  get_vars_code(d, vars, Float64)

# TODO vars can be vector of ... 
function get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol}, ::Type{stateeltype},) where {stateeltype}
  stmts = map(vars) do s
    ssymbl = QuoteNode(s)
    @match all(d[incident(d, s, :name), :type]) begin
      val && val .== :Constant => :($s = p.$s)
      val && val .== :Parameter => :($s = (p.$s)(t))
      val && val .== :Literal => :($s = $(parse(stateeltype, String(s))))
      _ => :($s = u.$s)
    end                 
    return quote
      $(stmts...)
    end
  end
end

# This is the setting of the du vector at the end of f
function set_tanvars_code(d::AbstractNamedDecapode)
  tanvars = [(d[e, [:src, :name]], d[e, [:tgt, :name]]) for e in incident(d, :∂ₜ, :op1)]
  stmts = map(tanvars) do (s, t)
    ssymb = QuoteNode(s)
    :(getproperty(du, $ssymb) .= $t)
  end
  return stmts
end

function compile_body!(visited::Dict, d::SummationDecapode, op::Int, ::Type{Val{:Op1}}) 
  if !visited[:Op1][op] && visited[:Var][d[op, :src]] && d[op, :op1] != DerivOp
    # skip the derivative edges
    operator = UnaryOperator(d[op, :op1])
    s, t = d[op, :src], d[op, :tgt]
    sname, tname = d[s, :name], d[t, :name]

    # equality = BinaryOperator(:(=)

    # TODO: Check to see if this is a DEC operator
    #
    if operator.is_optimizable && is_form(d, t)
      # TODO just set the is_broadcasted field to true 
      equality = promote_arithmetic_map(equality)
      operator = add_stub(in_place_stub, operator)

      push!(alloc_vectors, AllocVecCall(tname, d[t, :type], dimension, stateeltype))
    end
    call = UnaryCall(;operator=operator, equality=equality, input=sname, output=tname)

    # TODO we should instead have a single in-place function for setting an _operator_ to visited
    visited[:Op1][op], visited[:Var][t] = true, true
    
    call
  end
end

# TODO visited
function compile_body!(schedule::Vector{AbstractCall}, d::SummationDecapode, op::Int, ::Type{Val{:Op2}})

  # TODO this is a "table" in the ACSet
  r, arg1, arg2 = d[op, :res], d[op, :proj1], d[op, :proj2]
  # TODO these are foreign tables. can we use MLStyle to query these table effectively?
  a1name, a2name, rname = d[arg1, :name], d[arg2, :name], d[r, :name]

  equality = :(=)
  #operator = operator ∈ [:(*), :(-)] ? promote_arithmetic_map(d[op, :op2]) : d[op, :op2]
    operator = to_op(d, op, :op2)

  visited[:Op2][op], visited[:Var][t] = true, true
  
  # TODO unfinished
  call = BinaryCall(operator,
                    :(=),

  # the point of this block is to loop over d[:res] ∩ forms
  any_inferred = any(is_infer.(d, [arg1, arg2]))
  if is_form(d, r)
    @match operator begin
      operator.optype == :A || (operator.optype == :M && if any_inferred end) => begin
        operator.is_broadcasted = true
        # TODO: Repeated in both blocks:
        equality.is_broadcasted = true
        push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype))
      end
      op && if op.is_optimizable end => begin
        operator = add_stub(in_place_stub, op)
        # push!(operator.stub, in_place_stub)
        # TODO: Repeated in both blocks:
        equality.is_broadcasted = true
        push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype))
      end
      _ => nothing
    end
  end

  BinaryCall(operator, equality, a1name, a2name, rname)
  # TODO if we use a do-notation we can just return the call and schedule would accumulate them. we just need to modify the visited dictonary. that dictionary should just be internalized into the operators themselves;   operators should track if they've been visited.
end

function compile_body!(schedule::Vector{AbstractCall}, d::SummationDecapode, op::Int, ::Type{Val{:OpΣ}})
  r = d[op, :sum]
  summands = d[incident(d, op, :summation), :summand]
  argnames = d[summands, :name]
  rname = d[r, :name]

  operator = :(+)
  equality = :(=)

  # If result is a known form, broadcast addition
  # TODO: Also need to tell handler to prealloc
  if is_form(d, r)
    operator, equality = promote_arithmetic_map.(operator, equality)
    push!(alloc_vectors, AllocVecCall(rname, d[r, :type], dimension, stateeltype))
  end

  operator = :(.+)
  c = VarargsCall(operator, equality, argnames, rname)
  # visited_\Sigma \Var = true, true
  push!(schedule, c)
end

function compile(d::SummationDecapode, inputs::Vector, alloc_vectors::Vector{AllocVecCall}, optimizable_dec_operators::Set{Symbol};
  dimension = 2, stateeltype = Float64,)
  # Get the Vars of the inputs (probably state Vars).
  visited = merge([Dict(op => falses.(nparts(d, x))) for x ∈ [:Var, :Op1, :Op2, :Σ]]...)
  input_numbers = reduce(vcat, incident(d, inputs, :name))
  visited[:Var][input_numbers] .= true
  visited[:Var][incident(d, :Literal, :type)] .= true
  
  # FIXME: this is a quadratic implementation of topological_sort inlined in here.
  # TODO we should store in the terms whether they have been visited
  schedule = Vector{AbstractCall}()
  schedule = map(collect(Dict(:Op1 => :op1, :Op2 => :op2, :Σ => :Σ))) do (k,v)
    compile!.(visited, d, d[v], Val{k})
  end
  # e.g., d[:op1] = :Δ, :d, :∂ₜ, :∂q
  
  cache_exprs = map(alloc_vectors) do vec
    :($(vec.name) = (Decapodes.get_tmp($(Symbol(:__, vec.name)), u)))
  end

  eq_exprs = map(Expr, schedule)
  vcat(cache_exprs, eq_exprs)
end

# TODO: Add more specific types later for optimization
function resolve_types_compiler!(d::SummationDecapode)
  d[:type] = map(d[:type]) do x
    if (x == :Constant || x == :Parameter)
      return :infer
    end
    return x
  end
end

function replace_negation_with_multiply!(d::SummationDecapode, input_vars)
  found_negation = false
  rem_negations = []
  neg1var = 0

  # TODO unary operator can also say "is_negation"?
  # filter by operators which are negation
  map(enumerate(filter(op -> op.is_negation, d[:op1]))) do (i, op)
    ...
  end
  #
  for (i, op) in enumerate(d[:op1])
    if (op == :(-) || op == :neg)
      if (!found_negation)
        neg1var = add_part!(d, :Var, type = :Literal, name = Symbol("-1.0"))
        push!(input_vars, Symbol("-1.0"))
        found_negation = true
      end
      push!(rem_negations, i)
      add_part!(d, :Op2, proj1 = neg1var, proj2 = d[i, :src], res = d[i, :tgt], op2 = :.*)
    end

  end
  rem_parts!(d, :Op1, rem_negations)
end

function replace_names_compiler!(d::SummationDecapode)
  dec_op1, dec_op2 = Pair{Symbol,Any}[], Pair{Symbol,Symbol}[(:∧₀₀=>:.*)]
  replace_names!(d, dec_op1, dec_op2)
end

# TODO maybe the @match isn't necessary
function infer_overload_compiler!(d::SummationDecapode, dimension::Int)
  @match dimension begin
    1 => begin
      infer_types!(d, op1_inf_rules_1D, op2_inf_rules_1D)
      resolve_overloads!(d, op1_res_rules_1D, op2_res_rules_1D)
    end
    2 => begin
      infer_types!(d, op1_inf_rules_2D, op2_inf_rules_2D)
      resolve_overloads!(d, op1_res_rules_2D, op2_res_rules_2D)
    end
    _ => nothing
  end
end

# TODO I think out variables (or those being modified in place) should be first as arguments
function init_dec_matrices!(d::SummationDecapode, dec_matrices::Vector{Symbol}, optimizable_dec_operators::Set{Symbol})
  # TODO type issue with ∩ here but I think this is a virtuous suggestion
  for op ∈ union(d[:op1], d[:op2]) ∩ optimizable_dec_operators
    push!(dec_matrices, op)
  end
end

function link_contract_operators(d::SummationDecapode, con_dec_operators::Set{Symbol})

  contract_defs = quote end

  compute_to_name = Dict()
  curr_id = 1
  
  # wrong indexing but i feel we should be able to retrieve the table
  #
  map(p[:Op1] isa AbstractArray) do (id, name)
    computation = reverse!(map(x -> add_stub(InPlaceOperator(x)), name))
    compute_key = join(computation, " * ")

    computation_name = get(compute_to_name, compute_key, :Error)
    # if error

    d[id, :op1] = computation_name
  end

  contract_defs

    

  for (i, op) in enumerate(filter(op -> d[:op1] isa AbstractArray, d[:Op1]))


  for op1_id in parts(d, :Op1)
    op1_name = d[op1_id, :op1]
    if isa(op1_name, AbstractArray)
      computation = reverse!(map(x -> add_stub(InPlaceOperator(x)), op1_name))
      compute_key = join(computation, " * ")

      computation_name = get(compute_to_name, compute_key, :Error)
      if (computation_name == :Error)
        computation_name = add_stub(Symbol("GenSim-ConMat"), Symbol(curr_id))
        get!(compute_to_name, compute_key, computation_name)
        push!(con_dec_operators, computation_name)

        expr_line = Expr(
          Symbol("="),
          add_stub(computation_name, InPlaceStub),
          Expr(:call, :*, computation...),
        )
        push!(contract_defs.args, expr_line)

        expr_line = Expr(
          Symbol("="),
          computation_name,
          Expr(Symbol("->"), :x, Expr(:call, :*, add_inplace_stub(computation_name), :x)),
        )
        push!(contract_defs.args, expr_line)

        curr_id += 1
      end

      d[op1_id, :op1] = computation_name
    end
  end

  contract_defs
end

# TODO: Will want to eventually support contracted operations
function gensim(
  user_d::AbstractNamedDecapode,
  input_vars;
  dimension::Int = 2,
  stateeltype = Float64,
)
  # TODO: May want to move this after infer_types if we let users
  # set their own inference rules
  recognize_types(user_d)

  # Makes copy
  d′ = expand_operators(user_d)
  #d′ = average_rewrite(d′)

  replace_negation_with_multiply!(d′, input_vars)

  dec_matrices = Vector{Symbol}()
  alloc_vectors = Vector{AllocVecCall}()

  vars = get_vars_code(d′, input_vars, stateeltype)
  tars = set_tanvars_code(d′)

  # We need to run this after we grab the constants and parameters out
  resolve_types_compiler!(d′)
  infer_overload_compiler!(d′, dimension)

  # This should probably be followed by an expand_operators
  replace_names_compiler!(d′)
  open_operators!(d′, dimension = dimension)
  infer_overload_compiler!(d′, dimension)

  # This will generate all of the fundemental DEC operators present
  # TODO can we store "is_optimizable" in the operator structs?
  optimizable_dec_operators =
    Set([:⋆₀, :⋆₁, :⋆₂, :⋆₀⁻¹, :⋆₂⁻¹, :d₀, :d₁, :dual_d₀, :d̃₀, :dual_d₁, :d̃₁])
  extra_dec_operators = Set([:⋆₁⁻¹, :∧₀₁, :∧₁₀, :∧₁₁, :∧₀₂, :∧₂₀])

  init_dec_matrices!(
    d′,
    dec_matrices,
    union(optimizable_dec_operators, extra_dec_operators),
  )

  # This contracts matrices together into a single matrix
  contracted_dec_operators = Set{Symbol}()
  contract_operators!(d′, allowable_ops = optimizable_dec_operators)
  cont_defs = link_contract_operators(d′, contracted_dec_operators)

  union!(optimizable_dec_operators, contracted_dec_operators, extra_dec_operators)

  # Compilation of the simulation
  equations = compile(
    d′,
    input_vars,
    alloc_vectors,
    optimizable_dec_operators,
    dimension = dimension,
    stateeltype = stateeltype,
  )

  func_defs = compile_env(d′, dec_matrices, contracted_dec_operators)
  vect_defs = compile_var(alloc_vectors)

  quote
    function simulate(mesh, operators, hodge = GeometricHodge())
      $func_defs
      $cont_defs
      $vect_defs
      f(du, u, p, t) = begin
        $vars
        $(equations...)
        $(tars...)
      end
    end
  end
end

# TODO I think it is a moral issue that this function is homophonic with "gensym"
gensim(c::Collage; dimension::Int = 2) = gensim(collate(c); dimension = dimension)

"""    function gensim(d::AbstractNamedDecapode; dimension::Int=2)

Generate a simulation function from the given Decapode. The returned function can then be combined with a mesh and a function describing function mappings to return a simulator to be passed to `solve`.
"""
gensim(d::AbstractNamedDecapode; dimension::Int = 2, stateeltype = Float64) = gensim(
  d,
  vcat(collect(infer_state_names(d)), d[incident(d, :Literal, :type), :name]),
  dimension = dimension,
  stateeltype = stateeltype,
)

evalsim(d::AbstractNamedDecapode; dimension::Int = 2, stateeltype = Float64) =
  eval(gensim(d, dimension = dimension, stateeltype = stateeltype))
evalsim(d::AbstractNamedDecapode, input_vars; dimension::Int = 2, stateeltype = Float64) =
  eval(gensim(d, input_vars, dimension = dimension, stateeltype = stateeltype))

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

function closest_point(p1, p2, dims)
  p_res = collect(p2)
  for i = 1:length(dims)
    if dims[i] != Inf
      p = p1[i] - p2[i]
      f, n = modf(p / dims[i])
      p_res[i] += dims[i] * n
      if abs(f) > 0.5
        p_res[i] += sign(f) * dims[i]
      end
    end
  end
  Point3{Float64}(p_res...)
end

function flat_op(s::AbstractDeltaDualComplex2D, X::AbstractVector; dims = [Inf, Inf, Inf])
  # XXX: Creating this lookup table shouldn't be necessary. Of course, we could
  # index `tri_center` but that shouldn't be necessary either. Rather, we should
  # loop over incident triangles instead of the elementary duals, which just
  # happens to be inconvenient.
  tri_map = Dict{Int,Int}(triangle_center(s, t) => t for t in triangles(s))

  map(edges(s)) do e
    p = closest_point(point(s, tgt(s, e)), point(s, src(s, e)), dims)
    e_vec = (point(s, tgt(s, e)) - p) * sign(1, s, e)
    dual_edges = elementary_duals(1, s, e)
    dual_lengths = dual_volume(1, s, dual_edges)
    mapreduce(+, dual_edges, dual_lengths) do dual_e, dual_length
      X_vec = X[tri_map[s[dual_e, :D_∂v0]]]
      dual_length * dot(X_vec, e_vec)
    end / sum(dual_lengths)
  end
end
