begin
    using Decapodes
    using DiagrammaticEquations
    using CombinatorialSpaces
    using GeometryBasics
    using CUDA
    using LinearAlgebra
    using SparseArrays
    using BenchmarkTools
    Point2D = Point2{Float64}
    Point3D = Point3{Float64}
    CUDA.allowscalar(false)

  T = Float64
  
  # rect = loadmesh(Rectangle_30x10())
  rect = triangulated_grid(400, 400, 1, 1, Point3D)
  d_rect = EmbeddedDeltaDualComplex2D{Bool, T, Point3D}(rect)
  subdivide_duals!(d_rect, Circumcenter())
  
  R = zeros(T, ne(d_rect));
  R2 = zeros(T, ntriangles(d_rect));
  R3 = zeros(T, ntriangles(d_rect));

  A = rand(T, nv(d_rect));
  B = rand(T, ne(d_rect));
  B2 = rand(T, ne(d_rect));
  C = rand(T, ntriangles(d_rect))
  
  cuR = CuArray(R)
  cuR2 = CuArray(R2)
  cuR3 = CuArray(R3)
  cuA = CuArray(A)
  cuB = CuArray(B)
  cuB2 = CuArray(B2)
  cuC = CuArray(C)

  val_pack_01 = dec_p_wedge_product(Tuple{0,1}, d_rect)
  primal_vertices_01 = CuArray(val_pack_01[1])

  val_pack_02 = dec_p_wedge_product(Tuple{0,2}, d_rect)
  primal_vertices_02 = CuArray(val_pack_02[1])
  coeffs_02 = CuArray(val_pack_02[2])

  val_pack_11 = dec_p_wedge_product(Tuple{1,1}, d_rect)
  primal_vertices_11 = CuArray(val_pack_11[1])
  coeffs_11 = CuArray(val_pack_11[2])
end
  
  # Figure out a way to do atomic adds so we can split the addition across threads.y

  function dec_cu_c_wedge_product!(::Type{Tuple{0,1}}, wedge_terms, f, α, primal_vertices)
    num_threads = CUDA.max_block_size.x
    num_blocks = min(ceil(Int, length(wedge_terms) / num_threads), CUDA.max_grid_size.x)

    @cuda threads=num_threads blocks=num_blocks dec_cu_ker_c_wedge_product_01!(wedge_terms, f, α, primal_vertices)
  end

  function dec_cu_c_wedge_product!(::Type{Tuple{0,2}}, wedge_terms, f, α, primal_vertices, coeffs)
    num_threads = CUDA.max_block_size.x
    num_blocks = min(ceil(Int, length(wedge_terms) / num_threads), CUDA.max_grid_size.x)

    @cuda threads=num_threads blocks=num_blocks dec_cu_ker_c_wedge_product_02!(wedge_terms, f, α, primal_vertices, coeffs)
  end

  function dec_cu_c_wedge_product!(::Type{Tuple{1,1}}, wedge_terms, f, α, primal_vertices, coeffs)
    num_threads = CUDA.max_block_size.x
    num_blocks = min(ceil(Int, length(wedge_terms) / num_threads), CUDA.max_grid_size.x)

    @cuda threads=num_threads blocks=num_blocks dec_cu_ker_c_wedge_product_11!(wedge_terms, f, α, primal_vertices, coeffs)
  end

  function dec_cu_ker_c_wedge_product_01!(wedge_terms::CuDeviceArray{T}, f, α, primal_vertices) where T
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = gridDim().x * blockDim().x
    i = index
    @inbounds while i <= Int32(length(wedge_terms))
      wedge_terms[i] = T(0.5) * α[i] * (f[primal_vertices[i, Int32(1)]] + f[primal_vertices[i, Int32(2)]])
      i += stride
    end
    return nothing
  end

  function dec_cu_ker_c_wedge_product_02!(wedge_terms::CuDeviceArray{T}, f, α, pv, coeffs) where T
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = gridDim().x * blockDim().x
    i = index

    @inbounds while i <= Int32(length(wedge_terms))
      wedge_terms[i] = α[i] * (coeffs[Int32(1), i] * f[pv[Int32(1), i]] + coeffs[Int32(2), i] * f[pv[Int32(2), i]]
                                 + coeffs[Int32(3), i] * f[pv[Int32(3), i]] + coeffs[Int32(4), i] * f[pv[Int32(4), i]]
                                 + coeffs[Int32(5), i] * f[pv[Int32(5), i]] + coeffs[Int32(6), i] * f[pv[Int32(6), i]])
      i += stride
    end
    return nothing
  end

  function dec_cu_ker_c_wedge_product_11!(wedge_terms, α, β, e, coeffs)
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = gridDim().x * blockDim().x
    i = index

    @inbounds while i <= Int32(length(wedge_terms))
        e0, e1, e2 = e[Int32(1), i], e[Int32(2), i], e[Int32(3), i]
        ae0, ae1, ae2 = α[e0], α[e1], α[e2]
        be0, be1, be2 = β[e0], β[e1], β[e2]

        c1, c2, c3 = coeffs[Int32(1), i], coeffs[Int32(2), i], coeffs[Int32(3), i]

        wedge_terms[i] = (c1 * (ae2 * be1 - ae1 * be2)
                        + c2 * (ae2 * be0 - ae0 * be2)
                        + c3 * (ae1 * be0 - ae0 * be1))
        i += stride
    end

    return nothing
end

### Exterior Derivatives Here ###

function dec_p_derivbound(::Type{Val{0}}, sd::HasDeltaSet; transpose::Bool=false, negate::Bool=false)
  vec_size = 2 * ne(sd)

  # XXX: This is assuming that meshes don't have too many entries
  # TODO: This type should be settable by the user and default set to Int32
  I = Vector{Int32}(undef, vec_size)
  J = Vector{Int32}(undef, vec_size)

  V = Vector{Int8}(undef, vec_size)

  e_orient::Vector{Int8} = sd[:edge_orientation]
  for i in eachindex(e_orient)
      e_orient[i] = (e_orient[i] == 1 ? 1 : -1)
  end

  v0_list = @view sd[:∂v0]
  v1_list = @view sd[:∂v1]

  for i in edges(sd)
      j = 2 * i - 1

      I[j] = i
      I[j+1] = i

      J[j] = v0_list[i]
      J[j+1] = v1_list[i]

      sign_term = e_orient[i]

      V[j] = sign_term
      V[j+1] = -1 * sign_term
  end

  if (transpose)
      I, J = J, I
  end
  if (negate)
      V .= -1 .* V
  end

  (I, J, V)
end

function dec_differential_test!(::Type{Val{0}}, sd)
  _, J, V = dec_p_derivbound(Val{0}, sd)
  J_p = Array(reshape(J, 2, length(J) ÷ 2)')
  V_p = Array(reshape(V, 2, length(V) ÷ 2)')

  e_orient = @view sd[:edge_orientation]
  all_edge_orientations_same = all(e_orient .== e_orient[1])

  if(all_edge_orientations_same)
    if(e_orient[1] == -1)
      J_p .= J_p[:, [2,1]]
    end
    return (res, x) -> dec_c_differential_same!(res, x, J_p)
  else
    return (res, x) -> dec_c_differential_not_same!(res, x, J_p, V_p)
    return nothing
  end
end

function dec_c_differential_same!(res, f, indices)
  @inbounds for i in 1:length(res)
    ind1 = indices[i, 1]
    ind2 = indices[i, 2]
    res[i] = f[ind1] - f[ind2]
  end
end

function dec_c_differential_not_same!(res, f, indices, signs)
  @inbounds for i in 1:length(res)
    ind1 = indices[i, 1]
    ind2 = indices[i, 2]
    res[i] = signs[ind1] * (f[ind1] - f[ind2])
  end
end

function dec_differential_cu_test!(::Type{Val{0}}, sd)
  _, J, V = dec_p_derivbound(Val{0}, sd)
  J_p = CuArray(reshape(J, 2, length(J) ÷ 2)')
  V_p = CuArray(reshape(V, 2, length(V) ÷ 2)')

  e_orient = @view sd[:edge_orientation]
  all_edge_orientations_same = all(e_orient .== e_orient[1])

  if(all_edge_orientations_same)
    if(e_orient[1] == -1)
      J_p .= J_p[:, [2,1]]
    end
    return (res, x) -> dec_cu_c_differential_same!(res, x, J_p)
  else
    # return (res, x) -> dec_c_differential_not_same!(res, x, J_p, V_p)
  end
end

function dec_cu_c_differential_same!(res, x, J_p)
  num_threads = CUDA.max_block_size.x
  num_blocks = min(ceil(Int, length(res) / num_threads), CUDA.max_grid_size.x)

  @cuda threads=num_threads blocks=num_blocks dec_cu_ker_c_differential_same!(res, x, J_p)
end

function dec_cu_ker_c_differential_same!(res, f, indices)
  index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
  stride = gridDim().x * blockDim().x
  i = index

  @inbounds while i <= Int32(length(res))
    res[i] = f[indices[i, 1]] - f[indices[i, 2]]
    
    i += stride
  end

  return nothing
end



  #= 
  function dec_cu_atom_c_wedge_product_01!(wedge_terms::CuDeviceArray{T}, f, α, primal_vertices) where T
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = Int32(1024) * blockDim().x
    j = threadIdx().y
  
    i = index
    @inbounds while i <= length(wedge_terms)
      if(j == Int32(1))
        wedge_terms[i] = Int32(0)
      end
      sync_threads()
  
      temp = T(0.5) * α[i] * f[primal_vertices[i, j]]
      CUDA.atomic_add!(pointer(wedge_terms) + sizeof(T) * (i - Int32(1)), temp)
      sync_threads()
      
      i += stride
    end
    return nothing
  end
    
  function dec_cu_atom_share_c_wedge_product_01!(wedge_terms::CuDeviceArray{T}, f, α, primal_vertices) where T
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = Int32(1024) * blockDim().x
    j = threadIdx().y
  
    shared_arr = CUDA.CuDynamicSharedArray(T, blockDim().x)
  
    i = index
    @inbounds while i <= length(wedge_terms)
      if(threadIdx().y == Int32(1))
        shared_arr[threadIdx().x] = Int32(0)
      end
      sync_threads()
  
      temp = T(0.5) * α[i] * f[primal_vertices[i, j]]
      CUDA.atomic_add!(pointer(shared_arr) + sizeof(T) * (threadIdx().x - Int32(1)), temp)
      sync_threads()
      
      if(threadIdx().y == Int32(1))
        wedge_terms[i] = shared_arr[threadIdx().x]
      end
      sync_threads()
  
      i += stride
    end
    return nothing
  end
  
  function dec_cu_shuffle_c_wedge_product_01!(wedge_terms::CuDeviceArray{T}, f, α, primal_vertices) where T
    index = ((blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x)
    stride = Int32(512) * blockDim().x    
  
    i = index
    j = laneid() % Int32(2)
    @inbounds while i <= Int32(2 * length(wedge_terms))
        r_i = (i >> Int32(1)) + j
        val = T(0.5) * α[r_i] * f[primal_vertices[r_i, j + Int32(1)]]
        # @cuprintln("$(threadIdx().x), $(laneid()): ($(r_i),$(j+1)), $val")
        # mask = vote_ballot_sync(CUDA.FULL_MASK, (j == Int32(1)))
        val += shfl_up_sync(CUDA.FULL_MASK, val, Int32(1))
        # @cuprintln("$(threadIdx().x), $(laneid()): $val")

        if(j == Int32(0))
            wedge_terms[r_i] = val
        end

        i += stride
    end
    return nothing
  end
  =#