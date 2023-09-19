using CombinatorialSpaces
using LinearAlgebra
using Base.Iterators
using Catlab

function default_dec_matrix_generate(sd, my_symbol, hodge=GeometricHodge())
    op = @match my_symbol begin

        # Regular Hodge Stars
        :⋆₀ => dec_mat_hodge(0, sd, hodge)
        :⋆₁ => dec_mat_hodge(1, sd, hodge)
        :⋆₂ => dec_mat_hodge(2, sd, hodge)

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

        _ => error("Unmatched operator $my_symbol")
    end

    return op
end

function dec_mat_hodge(k, sd::HasDeltaSet, hodge)
    hodge = ⋆(k,sd,hodge=hodge)
    return (hodge, x-> hodge * x)
end

# TODO: Need to figure how to handle inverse hodge on
# DualForm1 in 2D due to it needing to take a matrix inverse
function dec_mat_inverse_hodge(k, sd::HasDeltaSet, hodge)
    invhodge = inv_hodge_star(k,sd,hodge)
    return (invhodge, x-> invhodge * x)
end

function dec_mat_differential(k, sd::HasDeltaSet)
    # diff = d(k,sd)
    diff = dec_p_differential(Val{k}, sd)
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
    
    op = @match my_symbol begin

        :plus => (+)
        :(-) || :neg => x-> -1 .* x
        :.* => (x,y) -> x .* y
        :./ => (x,y) -> x ./ y

        # Wedge products
        :∧₀₀ => dec_wedge_product(Tuple{0, 0}, sd)
        :∧₀₁ => dec_wedge_product(Tuple{0, 1}, sd)
        :∧₁₀ => dec_wedge_product(Tuple{1, 0}, sd)
        :∧₁₁ => dec_wedge_product(Tuple{1, 1}, sd)

        _ => default_dec_matrix_generate(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
end

# TODO: This relies on the assumption of a well ordering of the 
# the dual space simplices. If changed, use dec_p_wedge_product_zero
function dec_p_wedge_product_zero_one(sd)
    simples = simplices(1, sd)
    primal_vertices = map(x -> [sd[x, :∂v0], sd[x, :∂v1]], simples)
    return (primal_vertices, simples)
end

# TODO: This relies on the assumption of a well ordering of the 
# the dual space simplices. If changed, use dec_c_wedge_product_zero
# TODO: This assumes that the dual vertice on an edge is always the midpoint
function dec_c_wedge_product_zero_one(f, α, val_pack)
    primal_vertices, simples = val_pack

    wedge_terms = zeros(last(simples))
    @inbounds for i in simples
        wedge_terms[i] += f[primal_vertices[i][1]] + f[primal_vertices[i][2]]
    end

    return 0.5 .* wedge_terms .* α
end

function dec_p_wedge_product_zero(k, sd)

    # Gets a list of all of the 0 -> vertices, 1 -> edges, 2 -> triangles on mesh
    simples = simplices(k, sd)

    #These are a memory killers!!

    # For 1 -> edges, grabs the two dual edges that form the primal edge 
    # For 2 -> triangles, grabs all of the edges that radiate from the triangle center 
    subsimples = map(x -> subsimplices(k, sd, x), simples)

    # For 1 -> edges, gets the primal vertices of the dual edges 
    # For 2 -> triangles, gets primal vertices at the primal triangle corners
    primal_vertices = map(x -> primal_vertex(k, sd, x), subsimples)

    # Finding coeffs in wedge product is brutal on memory, around 345976 allocations for one map
    # vols = map(x -> volume(k,sd,x), simples)
    vols = CombinatorialSpaces.volume(k,sd,simples)
    dual_vols = map(y -> dual_volume(k,sd,y), subsimples)
    coeffs = dual_vols ./ vols
    return (primal_vertices, coeffs, simples)
end

# Remove any allocations for f_terms
function dec_c_wedge_product_zero(f, α, val_pack)
    primal_vertices, coeffs, simples = val_pack

    # TODO: May want to move this to be in the loop in case the coeffs width does change
    # Can use the below code in the preallocation to determine if we do have to recompute
    # the width at every step or if we can just compute it once.
    # all(map(x -> length(coeffs[x]), simples) .== length(coeffs[1]))
    width_iter = 1:length(coeffs[1])
    wedge_terms = zeros(last(simples))

    @inbounds for i in simples
        for j in width_iter
            wedge_terms[i] += coeffs[i][j] * f[primal_vertices[i][j]]
        end
    end
    
    return wedge_terms .* α
end

# This is adapted almost directly from the CombinatorialSpaces package
# Use this if some assumptions in the embedded delta sets changes
function dec_p_wedge_product_ones_safe(sd)
    simples = simplices(2, sd)

    coeffs = map(simples) do x
        dual_vs = vertex_center(sd, triangle_vertices(sd, x))
        dual_es = sort(incident(sd, triangle_center(sd, x), :D_∂v0),
                 by=e -> sd[e,:D_∂v1] .== dual_vs, rev=true)[1:3]
        map(dual_es) do e
            sum(dual_volume(2, sd, incident(sd, e, :D_∂e1)))
        end / CombinatorialSpaces.volume(2, sd, x)
    end
  
    e0 = map(x -> ∂(2,0,sd,x), simples)
    e1 = map(x -> ∂(2,1,sd,x), simples)
    e2 = map(x -> ∂(2,2,sd,x), simples)

    return (e0, e1, e2, coeffs, simples)
end

# TODO: This relies on a well established ordering for 
# the dual space simplices. If changed, use dec_p_wedge_product_ones_safe
function dec_p_wedge_product_ones(sd)
    simples = simplices(2, sd)

    coeffs = map(simples) do x
        dual_es = incident(sd, triangle_center(sd, x), :D_∂v0)[4:6]
        map(dual_es) do e
            sum(dual_volume(2, sd, incident(sd, e, :D_∂e1)))
        end / CombinatorialSpaces.volume(2, sd, x)
    end
  
    e0 = map(x -> ∂(2,0,sd,x), simples)
    e1 = map(x -> ∂(2,1,sd,x), simples)
    e2 = map(x -> ∂(2,2,sd,x), simples)

    return (e0, e1, e2, coeffs, simples)
end

function dec_c_wedge_product_ones(α, β, val_pack)
    e0, e1, e2, coeffs, simples = val_pack

    wedge_terms = zeros(last(simples))

    @inbounds for i in simples
        ae0, ae1, ae2 = α[e0[i]], α[e1[i]], α[e2[i]]
        be0, be1, be2 = β[e0[i]], β[e1[i]], β[e2[i]]

        wedge_terms[i] += (coeffs[i][1] * (ae2 * be1 - ae1 * be2) 
                         + coeffs[i][2] * (ae2 * be0 - ae0 * be2) 
                         + coeffs[i][3] * (ae1 * be0 - ae0 * be1))
    end

    return wedge_terms
end

function dec_wedge_product(::Type{Tuple{0,0}}, sd::HasDeltaSet)
    (f, g) -> f .* g
end

function dec_wedge_product(::Type{Tuple{1,0}}, sd::HasDeltaSet)
    val_pack = dec_p_wedge_product_zero_one(sd)
    (α, g) -> dec_c_wedge_product_zero_one(g, α, val_pack)
end

function dec_wedge_product(::Type{Tuple{0,1}}, sd::HasDeltaSet)
    val_pack = dec_p_wedge_product_zero_one(sd)
    (f, β) -> dec_c_wedge_product_zero_one(f, β, val_pack)
end

function dec_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product_zero(k, sd)
    (α, g) -> dec_c_wedge_product_zero(g, α, val_pack)
end

function dec_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product_zero(k, sd)
    (f, β) -> dec_c_wedge_product_zero(f, β, val_pack)
end

function dec_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D)
    val_pack = dec_p_wedge_product_ones(sd)
    (α, β) -> dec_c_wedge_product_ones(α, β,val_pack)
end

function default_dec_generate_1D(sd, my_symbol, hodge=GeometricHodge())
    
    op = @match my_symbol begin

        :L₀ => dec_lie_derivative_1D(Val{0}, sd, hodge)
        :L₁ => dec_lie_derivative_1D(Val{1}, sd, hodge)

        :i₁ => dec_interior_product_1D(Val{1}, sd, hodge)

        _ => default_dec_generate(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
end

function dec_interior_product_1D(::Type{Val{1}}, sd::HasDeltaSet, hodge)
    # Alpha is dual1, X is primal1

    # Takes alpha to primal0
    M_invhodge0 = inv_hodge_star(0,sd,hodge)
    invhodge0 = x -> M_invhodge0 * x

    # Alpha is primal0 and X is primal1, res is primal1
    wedge01 = dec_wedge_product(Tuple{0, 1}, sd)

    # Takes res to dual0
    M_hodge1 = ⋆(1,sd,hodge)
    hodge1 = x -> M_hodge1 * x

    (α, X) -> hodge1(wedge01(invhodge0(α), X))
end

function dec_lie_derivative_1D(::Type{Val{0}}, sd::HasDeltaSet, hodge)
    # Alpha is dual0, X is primal1

    # Takes alpha to dual1
    M_dual_d0 = dual_derivative(0,sd)
    dual_d0 = x-> M_dual_d0 * x

    #Alpha is dual1 and X is primal1, res is dual0 
    interior_d1_p1 = dec_interior_product_1D(Val{1}, sd, hodge)

    (X, α) -> interior_d1_p1(dual_d0(α), X)
end

function dec_lie_derivative_1D(::Type{Val{1}}, sd::HasDeltaSet, hodge)
    # Alpha is dual1, X is primal1

    #Alpha is dual1 and X is primal1, res is dual0 
    interior_d1_p1 = dec_interior_product_1D(Val{1}, sd, hodge)

    # Takes res to dual1
    M_dual_d0 = dual_derivative(0,sd)
    dual_d0 = x-> M_dual_d0 * x
    
    (X, α) -> dual_d0(interior_d1_p1(α, X))
end

function default_dec_generate_2D(sd, my_symbol, hodge=GeometricHodge())
    
    op = @match my_symbol begin

        :L₀ => dec_lie_derivative_2D(Val{0}, sd, hodge)
        :L₁ => dec_lie_derivative_2D(Val{1}, sd, hodge)
        :L₂ => dec_lie_derivative_2D(Val{2}, sd, hodge)

        :i₁ => dec_interior_product_2D(Val{1}, sd, hodge)
        :i₂ => dec_interior_product_2D(Val{2}, sd, hodge)

        _ => default_dec_generate(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
end

function dec_interior_product_2D(::Type{Val{1}}, sd::HasDeltaSet, hodge)
    # Alpha is dual1, X is primal1

    # Takes alpha to primal1
    if(hodge == GeometricHodge())
        M_hodge1 = -1 .* hodge_star(1,sd,hodge)
        invhodge1 = x -> (M_hodge1 \ x)
    else
        M_invhodge1 = inv_hodge_star(1,sd,hodge)
        invhodge1 = x -> M_invhodge1 * x
    end

    #Alpha term is primal1 and X are primal1, res is primal2
    wedge11 = dec_wedge_product(Tuple{1, 1}, sd)

    #Sends res to dual0
    M_hodge2 = -1 .* ⋆(2,sd,hodge)
    hodge2 = x -> M_hodge2 * x

    (α, X) -> hodge2(wedge11(invhodge1(α), X)) 
end

function dec_interior_product_2D(::Type{Val{2}}, sd::HasDeltaSet, hodge)
    # Alpha is dual2, X is primal1

    # Takes alpha to primal0
    M_invhodge0 = inv_hodge_star(0,sd,hodge)
    invhodge0 = x -> M_invhodge0 * x

    #Alpha term is primal0 and X are primal1, res is primal1
    wedge01 = dec_wedge_product(Tuple{0, 1}, sd)

    #Sends res to dual1
    M_hodge1 = ⋆(1,sd,hodge)
    hodge1 = x -> M_hodge1 * x

    (α, X) -> hodge1(wedge01(invhodge0(α), X))
end

function dec_lie_derivative_2D(::Type{Val{0}}, sd::HasDeltaSet, hodge)
    # Alpha is dual0, X is primal1

    # Takes alpha to dual1
    M_dual_d0 = dual_derivative(0,sd)
    dual_d0 = x-> M_dual_d0 * x

    # Alpha is dual1, X is primal1, res is dual0
    interior_d1_p1 = dec_interior_product_2D(Val{1}, sd, hodge)

    (X, α) -> interior_d1_p1(dual_d0(α), X)
end

function dec_lie_derivative_2D(::Type{Val{1}}, sd::HasDeltaSet, hodge)
    # Alpha is dual1, X is primal1

    # Takes alpha to dual2
    M_dual_d1 = dual_derivative(1,sd)
    dual_d1 = x-> M_dual_d1 * x

    # d(Alpha) is dual2, X is primal1, res_1 is dual1
    interior_d2_p1 = dec_interior_product_2D(Val{2}, sd, hodge)

    # Alpha is dual1, X is primal1, res_2 is dual0
    interior_d1_p1 = dec_interior_product_2D(Val{1}, sd, hodge)

    # Takes res_2 to dual1
    M_dual_d0 = dual_derivative(0,sd)
    dual_d0 = x-> M_dual_d0 * x

    (X, α) -> interior_d2_p1(dual_d1(α), X) + dual_d0(interior_d1_p1(α, X))
end

function dec_lie_derivative_2D(::Type{Val{2}}, sd::HasDeltaSet, hodge)
    # Alpha is dual2, X is primal1

    # Alpha is dual2, X is primal1, res is dual1
    interior_d2_p1 = dec_interior_product_2D(Val{2}, sd, hodge)

    #Sends res to dual1 then to dual2
    M_dual_d1 = dual_derivative(1,sd)
    dual_d1 = x -> M_dual_d1 * x
    
    (X, α) -> dual_d1(interior_d2_p1(α, X))
end

function dec_p_differential(::Type{Val{0}}, sd::HasDeltaSet)
    vec_size = 2 * ne(sd)

    I = zeros(Int64, vec_size)
    J = zeros(Int64, vec_size)
    V = zeros(Int64, vec_size)

    sign_term = 0

    for i in edges(sd)
        j = 2 * i - 1

        I[j] = i
        I[j + 1] = i

        sign_term = sign(1, sd, i)

        J[j] = sd[i, :∂v0]
        J[j + 1] = sd[i, :∂v1]

        V[j] = sign_term
        V[j + 1] = -sign_term
    end

    sparse(I, J, V)
end

function dec_p_differential(::Type{Val{1}}, sd::HasDeltaSet)
    vec_size = 3 * ntriangles(sd)

    I = zeros(Int64, vec_size)
    J = zeros(Int64, vec_size)
    V = zeros(Int64, vec_size)

    general_edge_sign = sign(1, sd, 1)
    no_recompute_signs = allequal(sd[:edge_orientation])

    for i in triangles(sd)
        j = 3 * i - 2

        I[j] = i
        I[j + 1] = i
        I[j + 2] = i

        tri_sign = sign(2, sd, i)

        J[j] = sd[i, :∂e0]
        J[j + 1] = sd[i, :∂e1]
        J[j + 2] = sd[i, :∂e2]

        edge_sign_0 = general_edge_sign
        edge_sign_1 = general_edge_sign
        edge_sign_2 = general_edge_sign
        
        if(!(no_recompute_signs))
            edge_sign_0 = sign(1, sd, J[j])
            edge_sign_1 = sign(1, sd, J[j + 1])
            edge_sign_2 = sign(1, sd, J[j + 2])
        end

        V[j] = edge_sign_0 * tri_sign
        V[j + 1] = -1 * edge_sign_1 * tri_sign
        V[j + 2] = edge_sign_2 * tri_sign
    end

    sparse(I, J, V)
end

function open_operators(d::SummationDecapode; dimension::Int = 2)
    e = deepcopy(d)
    open_operators!(e, dimension = dimension)
    return e  
end

function open_operators!(d::SummationDecapode; dimension::Int = 2)
    op2_remove_stack = Vector{Int}();
    for op2_idx in parts(d, :Op2)
        op2_proj1 = d[op2_idx, :proj1]
        op2_proj2 = d[op2_idx, :proj2]
        op2_res = d[op2_idx, :res]
        op2_name = d[op2_idx, :op2]

        remove_op2 = 0

        ## Make substitution for complex operator into components
        @match (op2_name, dimension) begin
            (:i₁ , 2) => begin remove_op2 = add_Inter_Prod_2D!(Val{1}, d, op2_proj1, op2_proj2, op2_res) end
            (:i₂ , 2) => begin remove_op2 = add_Inter_Prod_2D!(Val{2}, d, op2_proj1, op2_proj2, op2_res) end

            (:L₀, 2) => begin remove_op2 = add_Lie_2D!(Val{0}, d, op2_proj1, op2_proj2, op2_res) end
            (:L₁, 2) => begin remove_op2 = add_Lie_2D!(Val{1}, d, op2_proj1, op2_proj2, op2_res) end
            (:L₂, 2) => begin remove_op2 = add_Lie_2D!(Val{2}, d, op2_proj1, op2_proj2, op2_res) end
            _ => nothing
        end

        ## If sub was made, add the original operator to be removed
        (remove_op2 > 0) && push!(op2_remove_stack, op2_idx)
    end

    ## Remove all subbed operators
    rem_parts!(d, :Op2, op2_remove_stack)

    ## Infer types and resolves overloads for all newly subbed operators
    ## TODO: This can be removed by explicitly typing and overloading in the sub rules themselves
    if(dimension == 1)
        infer_types!(d, op1_inf_rules_1D, op2_inf_rules_1D)
        resolve_overloads!(d, op1_res_rules_1D, op2_res_rules_1D)
    elseif(dimension == 2)
        infer_types!(d, op1_inf_rules_2D, op2_inf_rules_2D)
        resolve_overloads!(d, op1_res_rules_2D, op2_res_rules_2D)
    end

    ## Add unique names for all newly added variables
    fill_names!(d, lead_symbol = Symbol("Gensim_Var_"));
end

function add_Inter_Prod_2D!(::Type{Val{1}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)

    ## Adds the hodge Dual1 to Primal1
    inv_hodge_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj1_Inter, tgt = inv_hodge_tgt, op1 = :⋆)

    ## Adds the wedge between Primal1 and Primal1
    wedge_res = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op2, proj1 = inv_hodge_tgt, proj2 = proj2_Inter, res = wedge_res, op2 = :∧)

    ## Adds the hodge Primal2 to Dual0
    hodge_res = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = wedge_res, tgt = hodge_res, op1 = :⋆)

    ## Outputs negated value
    add_part!(d, :Op1, src = hodge_res, tgt = res_Inter, op1 = :neg)
end

function add_Inter_Prod_2D!(::Type{Val{2}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)

    ## Adds the hodge Dual2 to Primal0
    inv_hodge_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj1_Inter, tgt = inv_hodge_tgt, op1 = :⋆)

    ## Adds the wedge between Primal0 and Primal1
    wedge_res = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op2, proj1 = inv_hodge_tgt, proj2 = proj2_Inter, res = wedge_res, op2 = :∧)

    ## Outputs result of hodge Primal1 to Dual1
    add_part!(d, :Op1, src = wedge_res, tgt = res_Inter, op1 = :⋆)
end

function add_Lie_2D!(::Type{Val{0}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Outputs result of dual derivative Dual0 to Dual1
    dual_d_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj2_Lie, tgt = dual_d_tgt, op1 = :d)

    ## Takes interior product of Primal1 and Dual1 to Dual0
    add_Inter_Prod_2D!(Val{1}, d, dual_d_tgt, proj1_Lie, res_Lie)
end

function add_Lie_2D!(::Type{Val{1}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Outputs result of dual derivative Dual1 to Dual2
    dual_d_1_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj2_Lie, tgt = dual_d_1_tgt, op1 = :d)

    ## Takes interior product of Primal1 and Dual2 to Dual1
    inter_product_2_res = add_part!(d, :Var, type = :infer, name = nothing)
    add_Inter_Prod_2D!(Val{2}, d, dual_d_1_tgt, proj1_Lie, inter_product_2_res)

    ## Takes interior product of Primal1 and Dual1 to Dual0
    inter_product_1_res = add_part!(d, :Var, type = :infer, name = nothing)
    add_Inter_Prod_2D!(Val{1}, d, proj2_Lie, proj1_Lie, inter_product_1_res)

    ## Outputs result of dual derivative Dual0 to Dual1
    dual_d_0_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = inter_product_1_res, tgt = dual_d_0_tgt, op1 = :d)

    ## Outputs sum of both dual_d_0 and inter_product_2
    add_part!(d, :Op2, proj1 = inter_product_2_res, proj2 = dual_d_0_tgt, res = res_Lie, op2 = :plus)
end

function add_Lie_2D!(::Type{Val{2}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Takes interior product of Primal1 and Dual2 to Dual1
    inter_product_2_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_Inter_Prod_2D!(Val{2}, d, proj1_Lie, proj2_Lie, inter_product_2_tgt)

    ## Outputs result of dual derivative Dual1 to Dual2
    add_part!(d, :Op1, src = inter_product_2_tgt, tgt = res_Lie, op1 = :d)
end
