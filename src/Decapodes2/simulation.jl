using MultiScaleArrays
using OrdinaryDiffEq
using GeometryBasics

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
    input
    output
end

Base.Expr(c::UnaryCall) = begin
    operator = c.operator
    if isa(c.operator, AbstractArray)
        operator = :(compose($operator))
    end
    return :($(c.output) = $operator($(c.input)))
end
                
struct BinaryCall <: AbstractCall 
    operator
    input1
    input2
    output
end

Base.Expr(c::BinaryCall) = begin
    operator = c.operator
    if isa(c.operator, AbstractArray)
        operator = :(compose($(c.operator)))
    end
    return :($(c.output) = $operator($(c.input1), $(c.input2)))
end

struct VarargsCall <: AbstractCall 
    operator
    inputs
    output
end

Base.Expr(c::VarargsCall) = begin
    operator = c.operator
    if isa(c.operator, AbstractArray)
        operator = :(compose($(c.operator)))
    end
    arglist = c.inputs
    return :($(c.output) = $operator($(arglist...)))
end

function infer_states(d::AbstractNamedDecapode)
    vars = map(parts(d, :Var)) do v
        if length(incident(d, v, :tgt)) == 0 &&
            length(incident(d, v, :res)) == 0
            # v isn't a derived value
            return v
        else
            return nothing
        end
    end
    return filter(!isnothing, vars)
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

function get_vars_code(d::AbstractNamedDecapode, vars::Vector{Symbol})
    stmts = map(vars) do s
        ssymbl = QuoteNode(s)
        if all(d[incident(d, s, :name) , :type] .== :Constant)
            :($s = p.$s)
        elseif all(d[incident(d, s, :name) , :type] .== :Parameter)
            :($s = (p.$s)(t))
        else
            :($s = findnode(u, $ssymbl).values)
        end
    end
    return quote $(stmts...) end
end


function set_tanvars_code(d::AbstractNamedDecapode, statevars::Vector{Symbol})
    tanvars = [(d[e, [:src,:name]], d[e, [:tgt,:name]]) for e in incident(d, :∂ₜ, :op1)]
    stmts = map(tanvars) do (s,t)
        ssymb = QuoteNode(s)
        :(findnode(du, $ssymb).values .= $t)
    end
    return quote $(stmts...) end
end

function compile_env(d::AbstractNamedDecapode)
  assumed_ops = Set([:+, :*, :-, :/])
  defs = quote end
  for op in d[:op1]
    if op == DerivOp
      continue
    end
    ops = QuoteNode(op)
    def = :($op = generate(mesh, $ops))
    push!(defs.args, def)
  end
  for op in d[:op2]
    if op in assumed_ops
      continue
    end
    ops = QuoteNode(op)
    def = :($op = generate(mesh, $ops))
    push!(defs.args, def)
  end
  return defs
end

gensim(d::AbstractNamedDecapode) = gensim(d, collect(infer_state_names(d)))

function gensim(d::AbstractNamedDecapode, input_vars)
  d′ = expand_operators(d)
  defs = compile_env(d′)
  rhs = compile(d′, input_vars)
  quote
    function simulate(mesh)
      $defs
      return $rhs
    end
  end
end

compile(d::AbstractNamedDecapode) = compile(d, infer_state_names(d))

function compile(d::NamedDecapode, inputs::Vector)
    input_numbers = incident(d, inputs, :name)
    visited = falses(nparts(d, :Var))
    visited[collect(flatten(input_numbers))] .= true
    consumed1 = falses(nparts(d, :Op1))
    consumed2 = falses(nparts(d, :Op2))
    # FIXME: this is a quadratic implementation of topological_sort inlined in here.
    op_order = []
    for iter in 1:(nparts(d, :Op1) + nparts(d,:Op2))
        for op in parts(d, :Op1)
            s = d[op, :src]
            if !consumed1[op] && visited[s]
                # skip the derivative edges
                operator = d[op, :op1]
                t = d[op, :tgt]
                if operator == DerivOp
                    continue
                end
                consumed1[op] = true
                visited[t] = true
                sname = d[s, :name]
                tname = d[t, :name]
                c = UnaryCall(operator, sname, tname)
                push!(op_order, c)
            end
        end

        for op in parts(d, :Op2)
            arg1 = d[op, :proj1]
            arg2 = d[op, :proj2]
            if !consumed2[op] && visited[arg1] && visited[arg2]
                r = d[op, :res]
                a1name = d[arg1, :name]
                a2name = d[arg2, :name]
                rname  = d[r, :name]
                operator = d[op, :op2]
                consumed2[op] = true
                visited[r] = true
                c = BinaryCall(operator, a1name, a2name, rname)
                push!(op_order, c)
            end
        end
    end
    assigns = map(Expr, op_order)
    ret = :(return)
    ret.args = d[d[:,:incl], :name]
    return quote f(du, u, p, t) = begin
        $(get_vars_code(d, inputs))
        $(assigns...)
        du .= 0.0
        $(set_tanvars_code(d, inputs))
    end; end
end

function compile(d::SummationDecapode, inputs::Vector)
    input_numbers = incident(d, inputs, :name)
    visited = falses(nparts(d, :Var))
    visited[collect(flatten(input_numbers))] .= true
    consumed1 = falses(nparts(d, :Op1))
    consumed2 = falses(nparts(d, :Op2))
    consumedΣ = falses(nparts(d, :Σ))
    # FIXME: this is a quadratic implementation of topological_sort inlined in here.
    op_order = []
    for iter in 1:(nparts(d, :Op1) + nparts(d,:Op2)) + nparts(d, :Σ)
        for op in parts(d, :Op1)
            s = d[op, :src]
            if !consumed1[op] && visited[s]
                # skip the derivative edges
                operator = d[op, :op1]
                t = d[op, :tgt]
                if operator == DerivOp
                    continue
                end
                consumed1[op] = true
                visited[t] = true
                sname = d[s, :name]
                tname = d[t, :name]
                c = UnaryCall(operator, sname, tname)
                push!(op_order, c)
            end
        end

        for op in parts(d, :Op2)
            arg1 = d[op, :proj1]
            arg2 = d[op, :proj2]
            if !consumed2[op] && visited[arg1] && visited[arg2]
                r = d[op, :res]
                a1name = d[arg1, :name]
                a2name = d[arg2, :name]
                rname  = d[r, :name]
                operator = d[op, :op2]
                consumed2[op] = true
                visited[r] = true
                c = BinaryCall(operator, a1name, a2name, rname)
                push!(op_order, c)
            end
        end

        for op in parts(d, :Σ)
            args = subpart(d, incident(d, op, :summation), :summand)
            if !consumedΣ[op] && all(visited[args])
                r = d[op, :sum]
                argnames = d[args, :name]
                rname  = d[r, :name]
                operator = :(+)
                consumedΣ[op] = true
                visited[r] = true
                c = VarargsCall(operator, argnames, rname)
                push!(op_order, c)
            end
        end
    end

    assigns = map(Expr, op_order)
    ret = :(return)
    ret.args = d[d[:,:incl], :name]
    return quote f(du, u, p, t) = begin
        $(get_vars_code(d, inputs))
        $(assigns...)
        du .= 0.0
        $(set_tanvars_code(d, inputs))
    end; end
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
