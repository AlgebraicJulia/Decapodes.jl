begin
    using Decapodes
    using DiagrammaticEquations
    using CombinatorialSpaces
    using GeometryBasics
    using MLStyle
    using ComponentArrays
    using CUDA
    using LinearAlgebra
    using SparseArrays
    Point2D = Point2{Float64}
    Point3D = Point3{Float64}
    CUDA.allowscalar(false)
  
  # rect = loadmesh(Rectangle_30x10())
  rect = triangulated_grid(400, 400, 1, 1, Point3D)
  d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(rect)
  subdivide_duals!(d_rect, Circumcenter())
  
  R = rand(ne(d_rect));
  A = rand(nv(d_rect));
  B = rand(ne(d_rect));
  
  cuR = CuArray(R)
  cuA = CuArray(A)
  cuB = CuArray(B)
  
  val_pack = dec_p_wedge_product(Tuple{0,1}, d_rect)
  primal_vertices = CuArray(val_pack[1])
end
  
  # Figure out a way to do atomic adds so we can split the addition across threads.y
  function dec_cu_c_wedge_product_01!(wedge_terms, f, α, primal_vertices)
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = Int32(512) * blockDim().x
    i = index
    @inbounds while i <= Int32(length(wedge_terms))
      wedge_terms[i] = 0.5 * α[i] * (f[primal_vertices[i, 1]] + f[primal_vertices[i, 2]])
      i += stride
    end
    return nothing
  end
  
  function dec_cu_atom_c_wedge_product_01!(wedge_terms::CuDeviceArray{T}, f, α, primal_vertices) where T
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
    stride = Int32(512) * blockDim().x
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
    stride = Int32(512) * blockDim().x
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
  