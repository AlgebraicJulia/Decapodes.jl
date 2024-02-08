using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using CUDA
using CairoMakie
using LinearAlgebra
using SparseArrays
Point2D = Point2{Float64}
Point3D = Point3{Float64}
CUDA.allowscalar(false)

rect = loadmesh(Rectangle_30x10())
# rect = triangulated_grid(100, 100, 1, 1, Point3D)
d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(rect)
subdivide_duals!(d_rect, Circumcenter())

Heat = @decapode begin
    U::Form0
    ∂ₜ(U) == 100 * Δ(U)
end

function simulate2(mesh, operators, hodge = GeometricHodge())
  begin      
      cu_A = CuArray(dec_inv_hodge_star(0, mesh, DiagonalHodge()))
      cu_B = CUDA.CUSPARSE.CuSparseMatrixCSC{Float64}(dec_differential(0, mesh))
      cu_C = CuArray(dec_hodge_star(1, mesh, DiagonalHodge()))
      cu_D = CUDA.CUSPARSE.CuSparseMatrixCSC{Float64}(dec_dual_derivative(1, mesh))
  end
  begin
      cu_mat_Δ = cu_A * cu_D * cu_C * cu_B
  end
  begin
      var2 = CuArray{Float64}(undef, nv(mesh))
      U̇ = CuArray{Float64}(undef, nv(mesh))
  end
  f(du, u, p, t) = begin
          begin
              U = u.U
              var"100" = 100.0
          end
          mul!(var2, cu_mat_Δ, U)
          U̇ .= var"100" .* var2
          getproperty(du, :U) .= U̇
      end
end

function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      x => error("Unmatched operator $my_symbol")
    end
    return op
  end
  
fₘ = simulate2(d_rect, generate)

U = CuArray(map(d_rect[:point]) do (x,_)
        return x
    end)

u₀ = ComponentArray(U=U)

constants_and_parameters = ()

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4))
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

CUDA.@time solve(prob, Tsit5());

begin
  frames = 100
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, rect, color=Vector{Float64}(soln(0).U), colormap=:jet, colorrange=extrema(soln(0).U))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "Heat_GPU.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    msh.color = Vector{Float64}(soln(t).U)
  end
end

# Figure out a way to do atomic adds so we can split the addition across threads.y
function dec_cu_c_wedge_product_01!(wedge_terms, f, α, primal_vertices)
  index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x   
  stride = Int32(4) * blockDim().x
  i = index
  @inbounds while i <= length(wedge_terms)
    wedge_terms[i] = 0.5f0 * α[i] * (f[primal_vertices[i, 1]] + f[primal_vertices[i, 2]])
    i += stride
  end
  return nothing
end