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
        :δ₁ => dec_mat_codifferential(1, sd, hodge)
        :δ₂ => dec_mat_codifferential(2, sd, hodge)

        # Laplace-de Rham
        :Δ₀ => dec_mat_laplace_de_rham(0, sd, hodge)
        :Δ₁ => dec_mat_laplace_de_rham(1, sd, hodge)
        :Δ₂ => dec_mat_laplace_de_rham(2, sd, hodge)

        # Wedge Products
        :∧₀₁ => dec_pair_wedge_product(Tuple{0, 1}, sd)
        :∧₁₀ => dec_pair_wedge_product(Tuple{1, 0}, sd)
        :∧₀₂ => dec_pair_wedge_product(Tuple{0, 2}, sd)
        :∧₂₀ => dec_pair_wedge_product(Tuple{2, 0}, sd)
        :∧₁₁ => dec_pair_wedge_product(Tuple{1, 1}, sd)

        _ => error("Unmatched operator $my_symbol")
    end

    return op
end

function dec_mat_hodge(k, sd::HasDeltaSet, hodge)
    hodge = dec_hodge_star(k, sd, hodge=hodge)
    return (hodge, x-> hodge * x)
end

# TODO: Need to figure how to handle inverse hodge on
# DualForm1 in 2D due to it needing to take a matrix inverse
function dec_mat_inverse_hodge(k::Int, sd::HasDeltaSet, hodge)
    invhodge = inv_hodge_star(k,sd,hodge)
    return (invhodge, x-> invhodge * x)
end

function dec_mat_differential(k::Int, sd::HasDeltaSet)
    diff = dec_differential(k, sd)
    return (diff, x-> diff * x)
end

function dec_mat_dual_differential(k::Int, sd::HasDeltaSet)
    dualdiff = dec_dual_derivative(k, sd)
    return (dualdiff, x-> dualdiff * x)
end

function dec_mat_codifferential(k::Int, sd::HasDeltaSet, hodge)
    codiff = δ(k, sd; hodge = hodge)
    return (codiff, x-> codiff * x)
end

function dec_mat_laplace_de_rham(k::Int, sd::HasDeltaSet, hodge)
    lpdr = dec_p_laplace_de_rham(k, sd, hodge)
    return (lpdr, x-> lpdr * x)
end

function dec_mat_laplace_beltrami(k::Int, sd::HasDeltaSet)
    lpbt = ∇²(k, sd)
    return (lpbt, x-> lpbt * x)
end

function dec_pair_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product(Tuple{0, k}, sd)
    ((y, α, g) -> dec_c_wedge_product!(Tuple{0, k}, y, g, α, val_pack), 
    (α, g) -> dec_c_wedge_product(Tuple{0, k}, g, α, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product(Tuple{0, k}, sd)
    ((y, f, β) -> dec_c_wedge_product!(Tuple{0, k}, y, f, β, val_pack), 
     (f, β) -> dec_c_wedge_product(Tuple{0, k}, f, β, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D)
    val_pack = dec_p_wedge_product(Tuple{1,1}, sd)
    ((y, α, β) -> dec_c_wedge_product!(Tuple{1,1}, y, α, β,val_pack), 
     (α, β) -> dec_c_wedge_product(Tuple{1,1}, α, β,val_pack))
end


function default_dec_generate(sd, my_symbol, hodge=GeometricHodge())
    
    op = @match my_symbol begin

        :plus => (+)
        :(-) || :neg => x-> -1 .* x
        :.* => (x,y) -> x .* y
        :./ => (x,y) -> x ./ y        

        # :⋆₁⁻¹ 

        _ => default_dec_matrix_generate(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
end

# TODO: This relies on the assumption of a well ordering of the 
# the dual space simplices. If changed, use dec_p_wedge_product_zero
function dec_p_wedge_product(::Type{Tuple{0, 1}}, sd)
    return (hcat(sd[:∂v0], sd[:∂v1]), simplices(1, sd))
end

# TODO: This relies on the assumption of a well ordering of the 
# the dual space simplices. If changed, use dec_c_wedge_product_zero
# TODO: This assumes that the dual vertice on an edge is always the midpoint
function dec_c_wedge_product!(::Type{Tuple{0, 1}}, wedge_terms, f, α, val_pack)
    primal_vertices, simples = val_pack

    # wedge_terms = Vector{Float64}(undef, last(simples))
    wedge_terms .= 0.5 .* α
    @inbounds for i in simples
        wedge_terms[i] *= (f[primal_vertices[i, 1]] + f[primal_vertices[i, 2]])
    end

    return wedge_terms
end

function dec_p_wedge_product(::Type{Tuple{0, 2}}, sd::HasDeltaSet)

    # Gets a list of all of the 0 -> vertices, 1 -> edges, 2 -> triangles on mesh
    simples = simplices(2, sd)

    #These are a memory killers!!

    # For 1 -> edges, grabs the two dual edges that form the primal edge 
    # For 2 -> triangles, grabs all of the edges that radiate from the triangle center 
    subsimples = map(x -> subsimplices(2, sd, x), simples)

    # For 1 -> edges, gets the primal vertices of the dual edges 
    # For 2 -> triangles, gets primal vertices at the primal triangle corners
    pv = primal_vertex(2, sd)
    primal_vertices = map(x -> pv[x], subsimples)

    # Finding coeffs in wedge product is brutal on memory, around 345976 allocations for one map
    # vols = map(x -> volume(k,sd,x), simples)
    vols = CombinatorialSpaces.volume(2,sd,simples)

    dv = sd[:dual_area]
    dual_vols = map(y -> dv[y], subsimples)

    coeffs = dual_vols ./ vols
    return (primal_vertices, coeffs, simples)
end

# Remove any allocations for f_terms
function dec_c_wedge_product!(::Type{Tuple{0, 2}}, wedge_terms, f, α, val_pack)
    primal_vertices, coeffs, simples = val_pack

    # TODO: May want to move this to be in the loop in case the coeffs width does change
    # Can use the below code in the preallocation to determine if we do have to recompute
    # the width at every step or if we can just compute it once.
    # all(map(x -> length(coeffs[x]), simples) .== length(coeffs[1]))
    @inbounds for i in simples
                pv = primal_vertices[i]
                coeff = coeffs[i]
                wedge_terms[i] = α[i] * (coeff[1] * f[pv[1]] + coeff[2] * f[pv[2]] + coeff[3] * f[pv[3]]
                                         + coeff[4] * f[pv[4]] + coeff[5] * f[pv[5]] + coeff[6] * f[pv[6]])
    end
    
    return wedge_terms
end

function dec_p_wedge_product(::Type{Tuple{0, k}}, sd) where k

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
function dec_c_wedge_product!(::Type{Tuple{0, k}}, wedge_terms, f, α, val_pack) where k
    primal_vertices, coeffs, simples = val_pack

    # TODO: May want to move this to be in the loop in case the coeffs width does change
    # Can use the below code in the preallocation to determine if we do have to recompute
    # the width at every step or if we can just compute it once.
    # all(map(x -> length(coeffs[x]), simples) .== length(coeffs[1]))
    width_iter = 1:length(coeffs[1])
    @inbounds for i in simples
                for j in width_iter
                    wedge_terms[i] += coeffs[i][j] * f[primal_vertices[i][j]]
                end
                wedge_terms[i] *= α[i]
    end
    
    return wedge_terms
end

# This is adapted almost directly from the CombinatorialSpaces package
# Use this if some assumptions in the embedded delta sets changes
function dec_p_wedge_product_safe(::Type{Tuple{1, 1}}, sd)
    simples = simplices(2, sd)

    coeffs = map(simples) do x
        dual_vs = vertex_center(sd, triangle_vertices(sd, x))
        dual_es = sort(incident(sd, triangle_center(sd, x), :D_∂v0),
                 by=e -> sd[e,:D_∂v1] .== dual_vs, rev=true)[1:3]
        map(dual_es) do e
            sum(dual_volume(2, sd, incident(sd, e, :D_∂e1)))
        end / CombinatorialSpaces.volume(2, sd, x)
    end
  
    e = Array{Int64}(undef, 3, last(simples))
    e[1, :] = ∂(2,0,sd)
    e[2, :] = ∂(2,1,sd)
    e[3, :] = ∂(2,2,sd)
    return (e, coeffs, simples)
end

# TODO: This relies on a well established ordering for 
# the dual space simplices. If changed, use dec_p_wedge_product_ones_safe
function dec_p_wedge_product(::Type{Tuple{1, 1}}, sd)
    simples = simplices(2, sd)

    coeffs::Vector{Vector{Float64}} = map(simples) do x
        dual_es = incident(sd, triangle_center(sd, x), :D_∂v0)[4:6]
        map(dual_es) do e
            sum(dual_volume(2, sd, incident(sd, e, :D_∂e1)))
        end / CombinatorialSpaces.volume(2, sd, x)
    end
  
    # e0 = ∂(2,0,sd)
    # e1 = ∂(2,1,sd)
    # e2 = ∂(2,2,sd)

    e = Array{Int64}(undef, 3, last(simples))
    e[1, :] = ∂(2,0,sd)
    e[2, :] = ∂(2,1,sd)
    e[3, :] = ∂(2,2,sd)
    return (e, coeffs, simples)
    # return (e0, e1, e2, coeffs, simples)

    # return(hcat(∂(2,0,sd), ∂(2,1,sd), ∂(2,2,sd)), coeffs, simples)
end

function dec_c_wedge_product!(::Type{Tuple{1, 1}}, wedge_terms, α, β, val_pack)
    # e0, e1, e2, coeffs, simples = val_pack
    e, coeffs, simples = val_pack

    for i in simples
        # ae0, ae1, ae2 = α[e0[i]], α[e1[i]], α[e2[i]]
        # be0, be1, be2 = β[e0[i]], β[e1[i]], β[e2[i]]

        ae0, ae1, ae2 = α[e[1, i]], α[e[2, i]], α[e[3, i]]
        be0, be1, be2 = β[e[1, i]], β[e[2, i]], β[e[3, i]]

        wedge_terms[i] += (coeffs[i][1] * (ae2 * be1 - ae1 * be2) 
                         + coeffs[i][2] * (ae2 * be0 - ae0 * be2) 
                         + coeffs[i][3] * (ae1 * be0 - ae0 * be1))
    end

    return wedge_terms
end

function dec_c_wedge_product(::Type{Tuple{m, n}}, f, α, val_pack) where {m, n}
    # The last item in the val_pack should always be the range of simplices
    wedge_terms = zeros(last(last(val_pack)))
    return dec_c_wedge_product!(Tuple{m, n}, wedge_terms, f, α, val_pack)
end

dec_wedge_product(m::Int, n::Int, sd::HasDeltaSet) = dec_wedge_product(Tuple{m,n}, sd::HasDeltaSet)

function dec_wedge_product(::Type{Tuple{0,0}}, sd::HasDeltaSet)
    (f, g) -> f .* g
end

function dec_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product(Tuple{0, k}, sd)
    (α, g) -> dec_c_wedge_product(Tuple{0, k}, g, α, val_pack)
end

function dec_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where k
    val_pack = dec_p_wedge_product(Tuple{0, k}, sd)
    (f, β) -> dec_c_wedge_product(Tuple{0, k}, f, β, val_pack)
end

function dec_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D)
    val_pack = dec_p_wedge_product(Tuple{1,1}, sd)
    (α, β) -> dec_c_wedge_product(Tuple{1,1}, α, β,val_pack)
end

function default_dec_generate_1D(sd, my_symbol, hodge=GeometricHodge())
    
    op = @match my_symbol begin

        _ => default_dec_generate(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
end

function default_dec_generate_2D(sd, my_symbol, hodge=GeometricHodge())
    
    op = @match my_symbol begin

        _ => default_dec_generate(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
end

# Boundary Operators
dec_boundary(n::Int, sd::HasDeltaSet) = sparse(dec_p_boundary(Val{n}, sd)...)

dec_p_boundary(::Type{Val{k}}, sd::HasDeltaSet; negate = false) where k = 
    dec_p_derivbound(Val{k - 1}, sd, transpose = true, negate = negate)

# Dual Derivative Operators
dec_dual_derivative(n::Int, sd::HasDeltaSet) = sparse(dec_p_dual_derivative(Val{n}, sd)...)

dec_p_dual_derivative(::Type{Val{0}}, sd::HasDeltaSet1D) = 
    dec_p_boundary(Val{1}, sd, negate = true)

dec_p_dual_derivative(::Type{Val{0}}, sd::HasDeltaSet2D) = 
    dec_p_boundary(Val{2}, sd)

dec_p_dual_derivative(::Type{Val{1}}, sd::HasDeltaSet2D) = 
    dec_p_boundary(Val{1}, sd, negate = true)

# Exterior Derivative Operators
dec_differential(n::Int, sd::HasDeltaSet) = sparse(dec_p_derivbound(Val{n}, sd)...)

function dec_p_derivbound(::Type{Val{0}}, sd::HasDeltaSet; transpose = false, negate = false)
    vec_size = 2 * ne(sd)

    I = Vector{Int64}(undef, vec_size)
    J = Vector{Int64}(undef, vec_size)
    V = Vector{Int64}(undef, vec_size)

    sign_term::Int = sign(1, sd, 1)
    e_orient = @view sd[:edge_orientation]
    recompute_signs::Bool = !(allequal(e_orient))

    # v0_list::Vector{Int64} = sd[:∂v0]
    # v1_list::Vector{Int64} = sd[:∂v1]

    v0_list = @view sd[:∂v0]
    v1_list = @view sd[:∂v1]

    for i in edges(sd)
        j = 2 * i - 1

        I[j] = i
        I[j + 1] = i

        J[j] = v0_list[i]
        J[j + 1] = v1_list[i]

        if(recompute_signs)
            sign_term = sign(1, sd, i)
        end

        V[j] = sign_term
        V[j + 1] = -1 * sign_term
    end
    
    if(transpose)
        I, J = J, I
    end
    if(negate)
        V .= -1 .* V
    end
    
    (I, J, V)
end

function dec_p_derivbound(::Type{Val{1}}, sd::HasDeltaSet; transpose = false, negate = false)
    vec_size = 3 * ntriangles(sd)

    I = Vector{Int64}(undef, vec_size)
    J = Vector{Int64}(undef, vec_size)
    V = Vector{Int64}(undef, vec_size)

    sign_term::Int = sign(1, sd, 1)
    e_orient = @view sd[:edge_orientation]
    recompute_signs::Bool = !(allequal(e_orient))
    # e0_list::Vector{Int64} = sd[:∂e0]
    # e1_list::Vector{Int64} = sd[:∂e1]
    # e2_list::Vector{Int64} = sd[:∂e2]

    e0_list = @view sd[:∂e0]
    e1_list = @view sd[:∂e1]
    e2_list = @view sd[:∂e2]

    tri_sign_list::Vector{Int64} = sign(2, sd)

    edge_sign_0::Int = sign_term
    edge_sign_1::Int = sign_term
    edge_sign_2::Int = sign_term

    for i in triangles(sd)
        j = 3 * i - 2

        I[j] = i
        I[j + 1] = i
        I[j + 2] = i

        J[j] = e0_list[i]
        J[j + 1] = e1_list[i]
        J[j + 2] = e2_list[i]

        if(recompute_signs)
            edge_sign_0 = sign(1, sd, J[j])
            edge_sign_1 = sign(1, sd, J[j + 1])
            edge_sign_2 = sign(1, sd, J[j + 2])
        end

        tri_sign = tri_sign_list[i]

        V[j] = edge_sign_0 * tri_sign
        V[j + 1] = -1 * edge_sign_1 * tri_sign
        V[j + 2] = edge_sign_2 * tri_sign
    end
    if(transpose)
        I, J = J, I
    end
    if(negate)
        V .= -1 .* V
    end

    (I, J, V)
end

# These are Diagonal Hodges 

# TODO: Check this Hodge with a 1D mesh
function dec_p_hodge_diag(::Type{Val{0}}, sd::AbstractDeltaDualComplex1D)
    num_v_sd = nv(sd)

    hodge_diag_0 = zeros(num_v_sd)

    v1_list = @view sd[:D_∂v1]
    dual_lengths = @view sd[:dual_length]
    
    for d_edge_idx in eachindex(v1_list)
        v1 = v1_list[d_edge_idx]
        if(1 <= v1 <= num_v_sd)
            hodge_diag_1[v1] += dual_lengths[d_edge_idx]
        end
    end
    return hodge_diag_0
end

# TODO: Check this Hodge with a 1D mesh
function dec_p_hodge_diag(::Type{Val{1}}, sd::AbstractDeltaDualComplex1D)
    return 1 ./ CombinatorialSpaces.volume(Val{1}, sd, edges(sd))
end

# TODO: This function could still use some work
#= function dec_p_hodge_diag(::Type{Val{0}}, sd::AbstractDeltaDualComplex2D)
    hodge_diag_0 = zeros(nv(sd))
    # centers = vertex_center(sd)
    dual_areas = sd[:dual_area]

    #This is the same as vcat(incident(sd, incident(sd, v, :D_∂v1), :D_∂e1)...)
    to_find = [:D_∂e1, :D_∂v1]
    for v in vertices(sd)
        # duals = incident(sd, centers[v], to_find)
        duals = incident(sd, v, to_find)
        # v_duals = incident(sd, v, :D_∂v1)
        # v_duals = v_duals[1:length(vduals) / 2 + 1]
        # duals = incident(sd, v_duals, :D_∂e1)
        for dual in duals
            hodge_diag_0[v] += dual_areas[dual]
        end
    end
    return hodge_diag_0
end =#

function dec_p_hodge_diag(::Type{Val{0}}, sd::AbstractDeltaDualComplex2D)
    hodge_diag_0 = zeros(nv(sd))

    dual_edges_1 = @view sd[:D_∂e1]
    dual_v_1 = @view sd[:D_∂v1]
    dual_areas = @view sd[:dual_area]

    for dual_tri in eachindex(dual_edges_1)
        v = dual_v_1[dual_edges_1[dual_tri]]
        hodge_diag_0[v] += dual_areas[dual_tri]
    end
    return hodge_diag_0
end

function dec_p_hodge_diag(::Type{Val{1}}, sd::AbstractDeltaDualComplex2D)
    num_v_sd = nv(sd)
    num_e_sd = ne(sd)

    hodge_diag_1 = zeros(num_e_sd)

    # v1_list::Vector{Int64} = sd[:D_∂v1]
    # dual_lengths::Vector{Float64} = sd[:dual_length]
    # lengths::Vector{Float64} = sd[:length]

    v1_list = @view sd[:D_∂v1]
    dual_lengths = @view sd[:dual_length]
    lengths = @view sd[:length]

    for d_edge_idx in eachindex(v1_list)
        v1_shift = v1_list[d_edge_idx] - num_v_sd
        if(1 <= v1_shift <= num_e_sd)
            hodge_diag_1[v1_shift] += dual_lengths[d_edge_idx] / lengths[v1_shift]
        end
    end
    return hodge_diag_1
end

function dec_p_hodge_diag(::Type{Val{2}}, sd::AbstractDeltaDualComplex2D)
    tri_areas = @view sd[:area]
    return 1 ./ tri_areas
end

dec_hodge_star(n::Int, sd::HasDeltaSet; hodge = GeometricHodge()) = dec_hodge_star(n, sd, hodge) 
dec_hodge_star(n::Int, sd::HasDeltaSet, ::DiagonalHodge) = dec_hodge_star(Val{n}, sd, DiagonalHodge())
dec_hodge_star(n::Int, sd::HasDeltaSet, ::GeometricHodge) = dec_hodge_star(Val{n}, sd, GeometricHodge())

dec_hodge_star(::Type{Val{k}}, sd::HasDeltaSet, ::DiagonalHodge) where k = 
    Diagonal(dec_p_hodge_diag(Val{k}, sd))

# These are Geometric Hodges 
# TODO: Still need implementation for Hodge 1 in 2D
dec_hodge_star(::Type{Val{0}}, sd::AbstractDeltaDualComplex1D, ::GeometricHodge) = 
    dec_hodge_star(Val{0}, sd, DiagonalHodge())

dec_hodge_star(::Type{Val{1}}, sd::AbstractDeltaDualComplex1D, ::GeometricHodge) = 
    dec_hodge_star(Val{1}, sd, DiagonalHodge())

dec_hodge_star(::Type{Val{0}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge) = 
    dec_hodge_star(Val{0}, sd, DiagonalHodge())

dec_hodge_star(::Type{Val{2}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge) = 
    dec_hodge_star(Val{2}, sd, DiagonalHodge())

function crossdot(a, b)
    x, y, z = 1, 2, 3
    c_x = a[y] * b[z] - a[z] * b[y]
    c_y = a[z] * b[x] - a[x] * b[z]
    c_z = a[x] * b[y] - a[y] * b[x]

    flipbit = (c_z == 0 ? 1.0 : sign(c_z))
    c_norm = sqrt(c_x^2 + c_y^2 + c_z^2)
    return c_norm * flipbit
end

function dec_hodge_star(::Type{Val{1}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge)

    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    # sizehint!(I, ne(sd))
    # sizehint!(J, ne(sd))
    # sizehint!(V, ne(sd))

    rel_orient::Float64 = 0.0
    # tri_edges_0 = sd[:∂e0]
    # tri_edges_1 = sd[:∂e1]
    # tri_edges_2 = sd[:∂e2]

    edge_centers = @view sd[:edge_center]
    tri_centers = @view sd[:tri_center]

    points::Vector{Point3{Float64}} = sd[:point]
    dual_points::Vector{Point3{Float64}} = sd[:dual_point]

    # points = sd[:point]
    # dual_points = sd[:dual_point]

    tgts = @view sd[:∂v0]
    srcs = @view sd[:∂v1]

    tri_signs::Vector{Int64} = sign(2, sd)

    for t in triangles(sd)
      e = reverse(triangle_edges(sd, t))
      # ev = point(sd, tgt(sd, e)) .- point(sd, src(sd,e))
      ev = points[tgts[e]] .- points[srcs[e]]

      # tc = dual_point(sd, triangle_center(sd, t))
      tc = dual_points[tri_centers[t]]
      #= dv = map(enumerate(dual_point(sd, edge_center(sd, e)))) do (i,v)
        (tc - v) * (i == 2 ? -1 : 1)
      end =#
  
      dv = map(enumerate(dual_points[edge_centers[e]])) do (i,v)
        (tc - v) * (i == 2 ? -1 : 1)
      end

      diag_dot = map(1:3) do i
        dot(ev[i], dv[i]) / dot(ev[i], ev[i])
      end
  
      # This relative orientation needs to be redefined for each triangle in the
      # case that the mesh has multiple independent connected components
      rel_orient = 0.0
      for i in 1:3
        diag_cross = tri_signs[t] * crossdot(ev[i], dv[i]) /
                        dot(ev[i], ev[i])
        if diag_cross != 0.0
          # Decide the orientation of the mesh relative to z-axis (see crossdot)
          # For optimization, this could be moved out of this loop
          if rel_orient == 0.0
            rel_orient = sign(diag_cross)
          end
  
          push!(I, e[i])
          push!(J, e[i])
          push!(V, diag_cross * rel_orient)
        end
      end
  
      for p ∈ ((1,2,3), (1,3,2), (2,1,3),
               (2,3,1), (3,1,2), (3,2,1))
        val = rel_orient * tri_signs[t] * diag_dot[p[1]] *
                dot(ev[p[1]], ev[p[3]]) / crossdot(ev[p[2]], ev[p[3]])
        if val != 0.0
          push!(I, e[p[1]])
          push!(J, e[p[2]])
          push!(V, val)
        end
      end
    end
    sparse(I,J,V)
  end
  

# These are Diagonal Inverse Hodges
function dec_inv_hodge(::Type{Val{k}}, sd::HasDeltaSet, ::DiagonalHodge) where k
    hdg = dec_p_hodge_diag(Val{k}, sd)
    mult_term = iseven(k*(ndims(sd)-k)) ? 1 : -1
    hdg .= (1 ./ hdg) .* mult_term
    return Diagonal(hdg)
end

# These are Geometric Inverse Hodges
dec_inv_hodge(::Type{Val{0}}, sd::AbstractDeltaDualComplex1D, ::GeometricHodge) = 
    dec_inv_hodge(Val{0}, sd, DiagonalHodge())

dec_inv_hodge(::Type{Val{1}}, sd::AbstractDeltaDualComplex1D, ::GeometricHodge) = 
    dec_inv_hodge(Val{1}, sd, DiagonalHodge())

dec_inv_hodge(::Type{Val{0}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge) = 
    dec_inv_hodge(Val{0}, sd, DiagonalHodge())

# TODO: Change this hodge to dec hodge when implemented
function dec_inv_hodge(::Type{Val{1}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge)
    hdg_lu = LinearAlgebra.factorize(dec_hodge_star(1, sd, GeometricHodge()))
    x -> hdg_lu \ x
end

dec_inv_hodge(::Type{Val{2}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge) = 
    dec_inv_hodge(Val{2}, sd, DiagonalHodge())

dec_p_laplace_de_rham(n::Int, sd::HasDeltaSet, hodge = GeometricHodge()) = dec_p_laplace_de_rham(Val{n}, sd, hodge)

dec_p_laplace_de_rham(::Type{Val{0}}, sd::HasDeltaSet, hodge = GeometricHodge()) = 
    return δ(1, sd; hodge = hodge) * dec_differential(0, sd)

dec_p_laplace_de_rham(::Type{Val{n}}, sd::HasDeltaSet, hodge = GeometricHodge()) where n = 
    return δ(n + 1, sd; hodge = hodge) * dec_differential(n, sd) + dec_differential(n - 1, sd) * δ(n, sd; hodge = hodge)

dec_p_laplace_de_rham(::Type{Val{1}}, sd::AbstractDeltaDualComplex1D, hodge = GeometricHodge()) = 
    return dec_differential(0, sd) * δ(1, sd; hodge = hodge)
    
dec_p_laplace_de_rham(::Type{Val{2}}, sd::AbstractDeltaDualComplex2D, hodge = GeometricHodge()) = 
    return dec_differential(1, sd) * δ(2, sd; hodge = hodge)

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
            (:i₁ , 1) => begin remove_op2 = add_Inter_Prod_1D!(Val{1}, d, op2_proj1, op2_proj2, op2_res) end

            (:i₁ , 2) => begin remove_op2 = add_Inter_Prod_2D!(Val{1}, d, op2_proj1, op2_proj2, op2_res) end
            (:i₂ , 2) => begin remove_op2 = add_Inter_Prod_2D!(Val{2}, d, op2_proj1, op2_proj2, op2_res) end

            (:L₀, 1) => begin remove_op2 = add_Lie_1D!(Val{0}, d, op2_proj1, op2_proj2, op2_res) end
            (:L₁, 1) => begin remove_op2 = add_Lie_1D!(Val{1}, d, op2_proj1, op2_proj2, op2_res) end

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

function add_Inter_Prod(d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
    ## Adds the hodge Dual to Primal
    inv_hodge_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj1_Inter, tgt = inv_hodge_tgt, op1 = :⋆)

    ## Adds the wedge between Primal and Primal
    wedge_res = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op2, proj1 = inv_hodge_tgt, proj2 = proj2_Inter, res = wedge_res, op2 = :∧)

    ## Adds the hodge Primal to Dual
    add_part!(d, :Op1, src = wedge_res, tgt = res_Inter, op1 = :⋆)
end

function add_Inter_Prod_1D!(::Type{Val{1}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
    add_Inter_Prod(d, proj1_Inter, proj2_Inter, res_Inter)
end

function add_Inter_Prod_2D!(::Type{Val{1}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
    ## Takes generic interior product
    pos_inter_prod = add_part!(d, :Var, type = :infer, name = nothing)
    add_Inter_Prod(d, proj1_Inter, proj2_Inter, pos_inter_prod)

    ## Outputs negated value
    add_part!(d, :Op1, src = pos_inter_prod, tgt = res_Inter, op1 = :neg)
end

function add_Inter_Prod_2D!(::Type{Val{2}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
    add_Inter_Prod(d, proj1_Inter, proj2_Inter, res_Inter)
end

function add_Lie_1D!(::Type{Val{0}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Outputs result of dual derivative Dual0 to Dual1
    dual_d_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj2_Lie, tgt = dual_d_tgt, op1 = :d)

    ## Takes interior product of Primal1 and Dual1 to Dual0
    add_Inter_Prod_1D!(Val{1}, d, dual_d_tgt, proj1_Lie, res_Lie)
end

function add_Lie_1D!(::Type{Val{1}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Takes interior product of Primal1 and Dual1 to Dual0
    inter_product_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_Inter_Prod_1D!(Val{1}, d, proj2_Lie, proj1_Lie, inter_product_tgt)

    ## Outputs result of dual derivative Dual0 to Dual1
    add_part!(d, :Op1, src = inter_product_tgt, tgt = res_Lie, op1 = :d)
end

function add_Lie_2D!(::Type{Val{0}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Outputs result of dual derivative Dual0 to Dual1
    dual_d_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_part!(d, :Op1, src = proj2_Lie, tgt = dual_d_tgt, op1 = :d)

    ## Takes interior product of Primal1 and Dual1 to Dual0
    add_Inter_Prod_1D!(Val{1}, d, dual_d_tgt, proj1_Lie, res_Lie)
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
    summation_tgt = add_part!(d, :Σ, sum = res_Lie)

    add_part!(d, :Summand, summand = inter_product_2_res, summation = summation_tgt)
    add_part!(d, :Summand, summand = dual_d_0_tgt, summation = summation_tgt)
end

function add_Lie_2D!(::Type{Val{2}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

    ## Takes interior product of Primal1 and Dual2 to Dual1
    inter_product_2_tgt = add_part!(d, :Var, type = :infer, name = nothing)
    add_Inter_Prod_2D!(Val{2}, d, proj2_Lie, proj1_Lie, inter_product_2_tgt)

    ## Outputs result of dual derivative Dual1 to Dual2
    add_part!(d, :Op1, src = inter_product_2_tgt, tgt = res_Lie, op1 = :d)
end
