using CombinatorialSpaces
using MultiScaleArrays
using OrdinaryDiffEq
using GeometryBasics
using LinearAlgebra
using Base.Iterators
using Decapodes
using Catlab.ACSetInterface
using MLStyle
import Catlab.Programs.GenerateJuliaPrograms: compile

struct VectorForm{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end

struct PhysicsState{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
    names::Vector{Symbol}
end

findname(u::PhysicsState, s::Symbol) = findfirst(isequal(s), u.names)
findnode(u::PhysicsState, s::Symbol) = u.nodes[findname(u, s)]

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
        Expr(:call, :mul!, c.output, operator, c.input)
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
end

# TODO: Check sizes on dual forms
Base.Expr(c::AllocVecCall) = begin
    resolved_form = @match c.form begin
        :Form0 => :V
        :Form1 => :E
        :Form2 => :Tri
        :DualForm0 => :Tri
        :DualForm1 => :E
        :DualForm2 => :V
        _ => return :AllocVecCall_Error
    end

    :($(c.name) = Vector{Float64}(undef, nparts(mesh, $(QuoteNode(resolved_form)))))
end

# TODO: Need to figure out how to deal with duals
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
    if(type == :Form0 || type == :Form1 || type == :Form2 || 
        type == :DualForm0 || type == :DualForm1 || type == :DualForm2)
        return true
    end
    return false
end

function is_form(d::SummationDecapode, var_name::Symbol)
    var_id = first(incident(d, var_name, :name))
    return is_form(d, var_id)
end

function is_literal(d::SummationDecapode, var_id::Int)
    return (d[var_id, :type] == :Literal)
end

function is_literal(d::SummationDecapode, var_name::Symbol)
    var_id = first(incident(d, var_name, :name))
    return is_literal(d, var_id)
end

function is_infer(d::SummationDecapode, var_id::Int)
    return (d[var_id, :type] == :infer)
end

function is_infer(d::SummationDecapode, var_name::Symbol)
    var_id = first(incident(d, var_name, :name))
    return is_infer(d, var_id)
end

function add_stub(stub_name::Symbol, var_name::Symbol)
    return Symbol("$(stub_name)_$(var_name)")
end


function infer_states(d::SummationDecapode)
    filter(parts(d, :Var)) do v
        length(incident(d, v, :tgt)) == 0 &&
        length(incident(d, v, :res)) == 0 &&
        length(incident(d, v, :sum)) == 0 &&
        d[v, :type] != :Literal
    end
end

infer_state_names(d) = d[infer_states(d), :name]

# This will be the function and matrix generation
function compile_env(d::AbstractNamedDecapode, dec_matrices::Vector{Symbol})
    assumed_ops = Set([:+, :*, :-, :/, :.+, :.*, :.-, :./])
    defined_ops = Set()

    defs = quote end

    for op in dec_matrices
        if(op in defined_ops)
            continue
        end

        quote_op = QuoteNode(op)
        mat_op = add_stub(:M, op)
        # Could turn this into a special call
        def = :(($mat_op, $op) = default_dec_matrix_generate(mesh, $quote_op, hodge))
        push!(defs.args, def)

        push!(defined_ops, op)
    end

    for op in d[:op1]
      if op == DerivOp
        continue
      end
      #= if typeof(op) <: AbstractArray
        for sub_op in op
          if(sub_op in defined_ops)
              continue
          end
  
          ops = QuoteNode(sub_op)
          def = :($sub_op = generate(mesh, $ops))
          push!(defs.args, def)
  
          push!(defined_ops, sub_op)
        end
        continue
      end =#
      if(op in defined_ops)
          continue
      end
  
      ops = QuoteNode(op)
      def = :($op = operators(mesh, $ops))
        
      push!(defs.args, def)
  
      push!(defined_ops, op)
    end
    for op in d[:op2]
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
    return quote $(map(Expr, alloc_vectors)...) end
end

# This is the block of parameter setting inside f
# TODO: Pass this an extra type parameter that sets the size of the Floats
function get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol})
    stmts = map(vars) do s
        ssymbl = QuoteNode(s)
        if all(d[incident(d, s, :name) , :type] .== :Constant)
        #if only(d[incident(d, s, :name) , :type]) == :Constant
            :($s = p.$s)
        elseif all(d[incident(d, s, :name) , :type] .== :Parameter)
            :($s = (p.$s)(t))
        elseif all(d[incident(d, s, :name) , :type] .== :Literal)
            # TODO: Fix this. We assume that all literals are Float64s.
            :($s = $(parse(Float64, String(s))))
        else
            # TODO: If names are not unique, then the type is assumed to be a
            # form for all of the vars sharing a same name.
            :($s = findnode(u, $ssymbl).values)
        end
    end
    return quote $(stmts...) end
end

# This is the setting of the du vector at the end of f
function set_tanvars_code(d::AbstractNamedDecapode)
    tanvars = [(d[e, [:src,:name]], d[e, [:tgt,:name]]) for e in incident(d, :∂ₜ, :op1)]
    stmts = map(tanvars) do (s,t)
        ssymb = QuoteNode(s)
        :(findnode(du, $ssymb).values .= $t)
    end
    return stmts
end

function compile(d::SummationDecapode, inputs::Vector, dec_matrices::Vector{Symbol}, alloc_vectors::Vector{AllocVecCall})
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

    promote_arithmetic_map = Dict(:(+) => :.+, :(-) => :.-, :(*) => :.*, :(/) => :./, :(=) => :.=,
                                  :.+ => :.+, :.- => :.-, :.* => :.*, :./ => :./, :.= => :.=)

    optimizable_dec_operators = Set([:⋆₀, :⋆₁, :⋆₂, :⋆₀⁻¹, :⋆₁⁻¹, 
                                    :d₀, :d₁, :dual_d₀, :d̃₀, :dual_d₁, :d̃₁,
                                    :δ₀, :δ₁,
                                    :Δ₀, :Δ₁, :Δ₂])

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
                    push!(dec_matrices, operator)

                    if(is_form(d, t))
                        equality = promote_arithmetic_map[equality]
                        operator = add_stub(:M, operator)

                        push!(alloc_vectors, AllocVecCall(tname, d[t, :type]))
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
                if(is_form(d, r))
                    if(operator == :(+) || operator == :(-) || operator == :.+ || operator == :.-)
                        operator = promote_arithmetic_map[operator]
                        equality = promote_arithmetic_map[equality]
                        push!(alloc_vectors, AllocVecCall(rname, d[r, :type]))
                    
                    # TODO: Do we want to support the ability of a user to use the backslash operator?
                    elseif(operator == :(*) || operator == :(/) || operator == :.* || operator == :./)
                        # WARNING: This part may break if we add more compiler types that have different 
                        # operations for basic and broadcast modes, e.g. matrix multiplication vs broadcast
                        if(!is_infer(d, arg1) && !is_infer(d, arg2))
                            operator = promote_arithmetic_map[operator]
                            equality = promote_arithmetic_map[equality]
                            push!(alloc_vectors, AllocVecCall(rname, d[r, :type]))
                        end
                    end

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
                    push!(alloc_vectors, AllocVecCall(rname, d[r, :type]))
                end

                visited_Σ[op] = true
                visited_Var[r] = true
                c = VarargsCall(operator, equality, argnames, rname)
                push!(op_order, c)
            end
        end
    end

    return map(Expr, op_order)
end
  
# TODO: Add more specific types later for optimization
function resolve_types_compiler!(d::SummationDecapode)
    d[:type] = map(d[:type]) do x
        if(x == :Constant || x == :Parameter)
            return :infer
        end
        return x
    end
    d
end

# TODO: Will want to eventually support contracted operations
function gensim(user_d::AbstractNamedDecapode, input_vars)
    # TODO: May want to move this after infer_types if we let users
    # set their own inference rules
    recognize_types(user_d)

    # Makes copy
    d′ = expand_operators(user_d)
    #d′ = average_rewrite(d′)

    # Mutates
    infer_types!(d′)
    resolve_overloads!(d′)

    dec_matrices = Vector{Symbol}();
    alloc_vectors = Vector{AllocVecCall}();

    vars = get_vars_code(d′, input_vars)
    tars = set_tanvars_code(d′)
    
    # We need to run this after we grab the constants and parameters out
    resolve_types_compiler!(d′)

    # rhs = compile(d′, input_vars)
    equations = compile(d′, input_vars, dec_matrices, alloc_vectors)

    func_defs = compile_env(d′, dec_matrices)
    vect_defs = compile_var(alloc_vectors)

    quote
        function simulate(mesh, operators, hodge=GeometricHodge())
            $func_defs
            $vect_defs
            f(du, u, p, t) = begin
                $vars
                $(equations...)
                $(tars...)
            end;
        end
    end
end

gensim(d::AbstractNamedDecapode) = gensim(d,
    vcat(collect(infer_state_names(d)), d[incident(d, :Literal, :type), :name]))

evalsim(d::AbstractNamedDecapode) = eval(gensim(d))
evalsim(d::AbstractNamedDecapode, input_vars) = eval(gensim(d, input_vars))


function default_dec_matrix_generate(sd, my_symbol, hodge=GeometricHodge())
    op = @match my_symbol begin

        # Regular Hodge Stars
        :⋆₀ => dec_mat_hodge(0, sd, hodge)
        :⋆₁ => dec_mat_hodge(1, sd, hodge)
        :⋆₂ => dec_mat_hodge(1, sd, hodge)

        # Inverse Hodge Stars
        :⋆₀⁻¹ => dec_mat_inverse_hodge(0, sd, hodge)
        :⋆₁⁻¹ => dec_mat_inverse_hodge(1, sd, hodge)

        # Differentials
        :d₀ => dec_mat_differential(0, sd)
        :d₁ => dec_mat_differential(1, sd)

        # Dual Differentials
        :dual_d₀ || :d̃₀ => dec_mat_dual_differential(0, sd)
        :dual_d₁ || :d̃₁ => dec_mat_dual_differential(1, sd)

        # Codifferential
        # TODO: Why do we have a matrix type parameter which is unused?
        :δ₀ => dec_mat_codifferential(0, sd, hodge)
        :δ₁ => dec_mat_codifferential(1, sd, hodge)

        # Laplace-de Rham
        :Δ₀ => dec_mat_laplace_de_rham(0, sd)
        :Δ₁ => dec_mat_laplace_de_rham(1, sd)
        :Δ₂ => dec_mat_laplace_de_rham(2, sd)

        x => error("Unmatched operator $my_symbol")
    end

    return op
end

function dec_mat_hodge(k, sd::HasDeltaSet, hodge)
    hodge = ⋆(k,sd,hodge=hodge)
    return (hodge, x-> hodge * x)
end

function dec_mat_inverse_hodge(k, sd::HasDeltaSet, hodge)
    invhodge = inv_hodge_star(k,sd,hodge)
    return (invhodge, x-> invhodge * x)
end

function dec_mat_differential(k, sd::HasDeltaSet)
    diff = d(k,sd)
    return (diff, x-> diff * x)
end

function dec_mat_dual_differential(k, sd::HasDeltaSet)
    dualdiff = dual_derivative(k,sd)
    return (dualdiff, x-> dualdiff * x)
end

function dec_mat_codifferential(k, sd::HasDeltaSet, hodge)
    codiff = δ(k, sd, hodge, nothing)
    return (codiff, x-> codiff * x)
end

function dec_mat_laplace_de_rham(k, sd::HasDeltaSet)
    lpdr = Δ(k, sd)
    return (lpdr, x-> lpdr * x)
end

function dec_mat_laplace_beltrami(k, sd::HasDeltaSet)
    lpbt = ∇²(k, sd)
    return (lpbt, x-> lpbt * x)
end

function default_dec_generate(sd, my_symbol, hodge=GeometricHodge())

    # TODO: Need to change to cached version
    i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
    
    op = @match my_symbol begin

        :plus => (+)
        :(-) => x-> -x
        :neg => x-> -x
        :.* => (x,y) -> x .* y
        :./ => (x,y) -> x ./ y

        # Regular Hodge Stars
        :⋆₀ => dec_hodge(0, sd, hodge)
        :⋆₁ => dec_hodge(1, sd, hodge)
        :⋆₂ => dec_hodge(2, sd, hodge)

        # Inverse Hodge Stars
        :⋆₀⁻¹ => dec_inverse_hodge(0, sd, hodge)
        :⋆₁⁻¹ => dec_inverse_hodge(1, sd, hodge)

        # Differentials
        :d₀ => dec_differential(0, sd)
        :d₁ => dec_differential(1, sd)

        # Dual Differentials
        :dual_d₀ || :d̃₀ => dec_dual_differential(0, sd)
        :dual_d₁ || :d̃₁ => dec_dual_differential(1, sd)

        # Codifferential
        # TODO: Why do we have a matrix type parameter which is unused?
        :δ₀ => dec_codifferential(0, sd, hodge)
        :δ₁ => dec_codifferential(1, sd, hodge)

        # Laplace-de Rham
        :Δ₀ => dec_laplace_de_rham(0, sd)
        :Δ₁ => dec_laplace_de_rham(1, sd)
        :Δ₂ => dec_laplace_de_rham(2, sd)

        # Wedge products
        :∧₀₀ => (f, g) -> wedge_product(Tuple{0,0}, sd, x, y)
        :∧₀₁ => dec_wedge_product(Tuple{0, 1}, sd)
        :∧₁₀ => dec_wedge_product(Tuple{1, 0}, sd)

        # Lie Derivative 0
        :L₀ => dec_lie_derivative_zero(sd, hodge)

        x => error("Unmatched operator $my_symbol")
    end

    return (args...) ->  op(args...)
end

function dec_hodge(k, sd::HasDeltaSet, hodge)
    hodge = ⋆(k,sd,hodge=hodge)
    x-> hodge * x
end

function dec_inverse_hodge(k, sd::HasDeltaSet, hodge)
    invhodge = inv_hodge_star(k,sd,hodge)
    x-> invhodge * x
end

function dec_differential(k, sd::HasDeltaSet)
    diff = d(k,sd)
    x-> diff * x
end

function dec_dual_differential(k, sd::HasDeltaSet)
    dualdiff = dual_derivative(k,sd)
    x-> dualdiff * x
end

function dec_codifferential(k, sd::HasDeltaSet, hodge)
    codiff = δ(k, sd, hodge, nothing)
    x -> codiff * x
end

function dec_laplace_de_rham(k, sd::HasDeltaSet)
    lpdr = Δ(k, sd)
    x -> lpdr * x
end

function dec_laplace_beltrami(k, sd::HasDeltaSet)
    lpbt = ∇²(k, sd)
    x -> lpbt * x
end

function dec_lie_derivative_zero(sd::HasDeltaSet, hodge)
    tmphodge1 = dec_hodge(1, sd, hodge)
    tmpwedge10 = dec_wedge_product(Tuple{1, 0}, sd)
    (v, x)-> tmphodge1(tmpwedge10(v, x))
end

function dec_p_wedge_product_zero(k, sd)

    # Gets a list of all of the 0 -> vertices, 1 -> edges, 2 -> triangles on mesh
    simples = simplices(k, sd)

    #These are a memory killers!!

    # For 1 -> edges, grabs the two dual edges that form the primal edge 
    # For 2 -> triangles, grabs all of the edges that radiate from the triangle center 
    subsimples = map(x -> subsimplices(k, sd, x), simples)

    # For 1 -> edges, gets the primal vertices of the dual edges 
    primal_vertices = map(x -> primal_vertex(k, sd, x), subsimples)

    # Finding coeffs in wedge product is brutal on memory, around 345976 allocations for one map
    #vols = map(x -> volume(k,sd,x), simples)
    vols = CombinatorialSpaces.volume(k,sd,simples)
    dual_vols = map(y -> dual_volume(k,sd,y), subsimples)
    coeffs = dual_vols ./ vols
    return (primal_vertices, coeffs)
end

function dec_c_wedge_product_zero(k, f, α, val_pack)
    primal_vertices, coeffs = val_pack
    f_terms = map(x -> f[x], primal_vertices)

    lhs = dot.(coeffs, f_terms)
    return (lhs .*  α) ./ factorial(k)
end

function dec_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product_zero(k, sd)
    (α, g) -> dec_c_wedge_product_zero(k, g, α, val_pack)
end

function dec_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product_zero(k, sd)
    (f, β) -> dec_c_wedge_product_zero(k, f, β, val_pack)
end

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

"""
    function recursive_delete_parents!(d::SummationDecapode, to_delete::Vector{Int64})

Delete the given nodes and their parents in the decapode, recursively.
"""
function recursive_delete_parents(d::SummationDecapode, to_delete::Vector{Int64})
  e = SummationDecapode{Any, Any, Symbol}()
  copy_parts!(e, d, (:Var, :TVar, :Op1, :Op2, :Σ, :Summand))
  isempty(to_delete) || recursive_delete_parents!(e, to_delete)
  return e
end

function recursive_delete_parents!(d::SummationDecapode, to_delete::Vector{Int64})
  # TODO: We assume to_delete vars have no children. Is that okay?
  # TODO: Explicitly check that a Var in to_delete is not a summand.
  vars_to_remove = Vector{Int64}()
  s = Stack{Int64}()
  foreach(v -> push!(s, v), to_delete)
  while true
    curr = pop!(s)
    parents = reduce(vcat,
                     [d[incident(d, curr, :tgt), :src],
                      d[incident(d, curr, :res), :proj1],
                      d[incident(d, curr, :res), :proj2],
                      d[incident(d, curr, [:summation, :sum]), :summand]])

    # Remove the operations which have curr as the result.
    rem_parts!(d, :TVar,    incident(d, curr, :incl))
    rem_parts!(d, :Op1,     incident(d, curr, :tgt))
    rem_parts!(d, :Op2,     incident(d, curr, :proj1))
    rem_parts!(d, :Op2,     incident(d, curr, :proj2))
    rem_parts!(d, :Op2,     incident(d, curr, :res))
    # Note: Delete Summands before Σs.
    rem_parts!(d, :Summand, incident(d, curr, [:summation, :sum]))
    rem_parts!(d, :Σ,       incident(d, curr, :sum))

    # Do not remove parents which are used in some other computation. We rely
    # on the fact that a parent is guaranteed to point to curr.
    filter!(parents) do p
      # p must not be the src of any Op1.
      isempty(incident(d, p, :src)) &&
      # p must not be a proj of any Op2.
      isempty(incident(d, p, :proj1)) &&
      isempty(incident(d, p, :proj2)) &&
      # p must not be a summand of any summation.
      isempty(incident(d, p, :summand))
    end
    foreach(p -> push!(s, p), parents)

    push!(vars_to_remove, curr)

    isempty(s) && break
  end
  rem_parts!(d, :Var, sort!(unique!(vars_to_remove)))
end

function closest_point(p1, p2, dims)
    p_res = collect(p2)
    for i in 1:length(dims)
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

function flat_op(s::AbstractDeltaDualComplex2D, X::AbstractVector; dims=[Inf, Inf, Inf])
  # XXX: Creating this lookup table shouldn't be necessary. Of course, we could
  # index `tri_center` but that shouldn't be necessary either. Rather, we should
  # loop over incident triangles instead of the elementary duals, which just
  # happens to be inconvenient.
  tri_map = Dict{Int,Int}(triangle_center(s,t) => t for t in triangles(s))

  map(edges(s)) do e
    p = closest_point(point(s, tgt(s,e)), point(s, src(s,e)), dims)
    e_vec = (point(s, tgt(s,e)) - p) * sign(1,s,e)
    dual_edges = elementary_duals(1,s,e)
    dual_lengths = dual_volume(1, s, dual_edges)
    mapreduce(+, dual_edges, dual_lengths) do dual_e, dual_length
      X_vec = X[tri_map[s[dual_e, :D_∂v0]]]
      dual_length * dot(X_vec, e_vec)
    end / sum(dual_lengths)
  end
end
