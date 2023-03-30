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

# TODO: Need to figure out how to deal with duals
function get_form_number(d::SummationDecapode, var_id::Int)
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
end

function is_form(d::SummationDecapode, var_id::Int)
    return (get_form_number(d, var_id) != -1)
end

function is_form(d::SummationDecapode, var_name::Symbol)
    var_id = first(incident(d, var_name, :name))
    return (get_form_number(d, var_id) != -1)
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
    vars = map(parts(d, :Var)) do v
        if length(incident(d, v, :tgt)) == 0 &&
            length(incident(d, v, :res)) == 0 &&
            length(incident(d, v, :sum)) == 0
            # v isn't a derived value
            return v
        else
            return nothing
        end
    end
    return filter(!isnothing, vars)
end

infer_state_names(d) = d[infer_states(d), :name]

# This will be the function generation and the matrix, vector preallocation
function compile_env(d::AbstractNamedDecapode, dec_matrices::Vector{Symbol})
    assumed_ops = Set([:+, :*, :-, :/, :.+, :.*, :.-, :./])
    defined_ops = Set()
  
    defs = quote end
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

# This is the block of parameter setting inside f
# We need to pass this an extra type parameter that sets the type of the floats
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
    # return quote $(stmts...) end
    return stmts
end

function compile(d::SummationDecapode, inputs::Vector, dec_matrices::Vector{Symbol})
    # Get the Vars of the inputs (probably state Vars).
    visited_Var = falses(nparts(d, :Var))

    input_numbers = incident(d, inputs, :name)
    visited_Var[collect(flatten(input_numbers))] .= true
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

                # TODO: Check to see if this is a DEC operator
                if(operator in optimizable_dec_operators)
                    push!(dec_matrices, operator)

                    if(is_form(d, t))
                        equality = promote_arithmetic_map[equality]
                        operator = add_stub(:M, operator)
                    end
                end

                sname = d[s, :name]
                tname = d[t, :name]

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

                # TODO: Add call to preallocate result
                # TODO: Check to make sure that this logic never breaks
                if(is_form(d, r))
                    if(operator == :(+) || operator == :(-) || operator == :.+ || operator == :.-)
                        operator = promote_arithmetic_map[operator]
                        equality = promote_arithmetic_map[equality]

                    elseif(operator == :(*) || operator == :(/) || operator == :.* || operator == :./)
                        # WARNING: This part may break if we add more compiler types that have different 
                        # operations for basic and broadcast modes, e.g. matrix multiplication vs broadcast
                        if(!is_infer(d, arg1) && !is_infer(d, arg2))
                            operator = promote_arithmetic_map[operator]
                            equality = promote_arithmetic_map[equality]
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
                end

                visited_Σ[op] = true
                visited_Var[r] = true
                c = VarargsCall(operator, equality, argnames, rname)
                push!(op_order, c)
            end
        end
    end

    return assigns = map(Expr, op_order)
end

function recognize_types(d::AbstractNamedDecapode)
    unrecognized_types = setdiff(d[:type], [:Form0, :Form1, :Form2, :DualForm0,
                            :DualForm1, :DualForm2, :Literal, :Parameter,
                            :Constant, :infer])
    isempty(unrecognized_types) ||
      error("Types $unrecognized_types are not recognized.")
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
function gensim(d::AbstractNamedDecapode, input_vars)
    # TODO: May want to move this after infer_types if we let users
    # set their own inference rules
    recognize_types(d)

    # Makes copy
    d′ = expand_operators(d)
    #d′ = average_rewrite(d′)

    # Mutates
    infer_types!(d′)
    resolve_overloads!(d′)

    dec_matrices = Vector{Symbol}();

    vars = get_vars_code(d′, input_vars)
    tars = set_tanvars_code(d′)
    # We need to run this after we grab the constants and parameters out
    resolve_types_compiler!(d)

    # rhs = compile(d′, input_vars)
    equations = compile(d′, input_vars, dec_matrices)
    println(dec_matrices)
    defs = compile_env(d′, dec_matrices)
    quote
        function simulate(mesh, operators)
            $defs
            f(du, u, p, t) = begin
                $vars
                $(equations...)
                # du .= 0.0
                $(tars...)
            end;
        end
    end
end

gensim(d::AbstractNamedDecapode) = gensim(d,
    vcat(collect(infer_state_names(d)), d[incident(d, :Literal, :type), :name]))

function compiler_dec_matrix_generate(sd, my_symbol, hodge=GeometricHodge())
    
    op = @match my_symbol begin

        # Regular Hodge Stars
        :⋆₀ => dec_hodge(0, sd, hodge)
        :⋆₁ => dec_hodge(1, sd, hodge)
        :⋆₂ => dec_hodge(1, sd, hodge)

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

        x => error("Unmatched operator $my_symbol")
    end

    return (args...) ->  op(args...)
end

function dec_hodge(k, sd::HasDeltaSet, hodge)
    hodge = ⋆(k,sd,hodge=hodge)
    return (hodge, x-> hodge * x)
end

function dec_inverse_hodge(k, sd::HasDeltaSet, hodge)
    invhodge = inv_hodge_star(k,sd,hodge)
    return (invhodge, x-> invhodge * x)
end

function dec_differential(k, sd::HasDeltaSet)
    diff = d(k,sd)
    return (diff, x-> diff * x)
end

function dec_dual_differential(k, sd::HasDeltaSet)
    dualdiff = dual_derivative(k,sd)
    return (dualdiff, x-> dualdiff * x)
end

function dec_codifferential(k, sd::HasDeltaSet, hodge)
    codiff = δ(k, sd, hodge, nothing)
    return (codiff, x-> codiff * x)
end

function dec_laplace_de_rham(k, sd::HasDeltaSet)
    lpdr = Δ(k, sd)
    return (lpdr, x-> lpdr * x)
end

function dec_laplace_beltrami(k, sd::HasDeltaSet)
    lpbt = ∇²(k, sd)
    return (lpbt, x-> lpbt * x)
end