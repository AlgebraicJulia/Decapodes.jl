using MultiScaleArrays
using OrdinaryDiffEq
using GeometryBasics
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
    input
    output
end

Base.Expr(c::UnaryCall) = begin
    operator = c.operator
    if isa(c.operator, AbstractArray)
        operator = Expr(:call, :∘, reverse(operator)...)
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
  defined_ops = Set()

  defs = quote end
  for op in d[:op1]
    if op == DerivOp
      continue
    end
    if typeof(op) <: AbstractArray
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
    end
    if(op in defined_ops)
        continue
    end

    ops = QuoteNode(op)
    def = :($op = generate(mesh, $ops))
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

gensim(d::AbstractNamedDecapode) = gensim(d, collect(infer_state_names(d)))

function gensim(d::AbstractNamedDecapode, input_vars)
  #d′ = expand_operators(d)
  d′ = d
  #d′ = average_rewrite(d′)
  defs = compile_env(d′)
  rhs = compile(d′, input_vars)
  quote
    function simulate(mesh, operators)
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

function default_dec_generate(sd, my_symbol, hodge=GeometricHodge())

    # TODO: Need to change to cahced version
    i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
    
    op = @match my_symbol begin

        :plus => (+)
        :(-) => x-> -x
        :.* => (x,y) -> x .* y
        :./ => (x,y) -> x ./ y

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
