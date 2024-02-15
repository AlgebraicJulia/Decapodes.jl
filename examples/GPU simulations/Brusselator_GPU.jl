using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using CairoMakie
using GeometryBasics: Point2, Point3
using CUDA
using CUDA.CUSPARSE
Point2D = Point2{Float64}
Point3D = Point3{Float64}

Brusselator = @decapode begin
  # Values living on vertices.
  (U, V)::Form0{X} # State variables.
  (U2V)::Form0{X} # Named intermediate variables.
  (U̇, V̇)::Form0{X} # Tangent variables.
  # Scalars.
  (α)::Constant{X}
  F::Parameter{X}
  # A named intermediate variable.
  U2V == (U .* U) .* V
  # Specify how to compute the tangent variables.
  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(U))
  # Associate tangent variables with a state variable.
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

# TODO: Create square domain of approximately 32x32 vertices.
# s = loadmesh(Rectangle_30x10())
# Visualize the mesh.
s = triangulated_grid(200, 200, 1, 1, Point3D)
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0])
s[:point] = map(x -> scaling_mat*x, s[:point])
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

U = CuVector(map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end)

V = CuVector(map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end)

F₁ = CuVector(map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end)

F₂ = CUDA.zeros(nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₂ : F₁)

function sim(mesh, operators, hodge = GeometricHodge())
  begin
      (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_matrix_generate(mesh, :dual_d₁, hodge)
      (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_matrix_generate(mesh, :⋆₀⁻¹, hodge)
      (var"GenSim-M_d₀", d₀) = default_dec_matrix_generate(mesh, :d₀, hodge)
      (var"GenSim-M_⋆₁", ⋆₁) = default_dec_matrix_generate(mesh, :⋆₁, hodge)
  end
  begin
      var"GenSim-M_GenSim-ConMat_1" = CuSparseMatrixCSC(var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁" * var"GenSim-M_d₀")
  end
  begin
      sum_2 = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•6" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•2" = CuVector{Float64}(undef, nparts(mesh, :V))
      U2V = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•4" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•5" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•8" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•7" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•1" = CuVector{Float64}(undef, nparts(mesh, :V))
      V̇ = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•3" = CuVector{Float64}(undef, nparts(mesh, :V))
      U̇ = CuVector{Float64}(undef, nparts(mesh, :V))
  end
  f(du, u, p, t) = begin
          begin
              U = u.U
              V = u.V
              α = p.α
              F = p.F(t)
              var"1" = 1.0
              var"4.4" = 4.4
              var"3.4" = 3.4
          end
          mul!(sum_2, var"GenSim-M_GenSim-ConMat_1", U)
          mul!(var"•6", var"GenSim-M_GenSim-ConMat_1", U)
          var"•2" .= U .* U
          U2V .= var"•2" .* V
          var"•4" .= var"4.4" .* U
          var"•5" .= α .* var"•6"
          var"•8" .= var"3.4" .* U
          var"•7" .= var"•8" .- U2V
          var"•1" .= α .* sum_2
          sum_1 = (.+)(var"1", U2V)
          V̇ .= (.+)(var"•7", var"•1")
          var"•3" .= sum_1 .- var"•4"
          U̇ .= (.+)(var"•3", var"•5", F)
          getproperty(du, :U) .= U̇
          getproperty(du, :V) .= V̇
      end
end

fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=V)

tₑ = 25

begin
  @info("Precompiling Solver")
  prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
  soln = solve(prob, Tsit5())
  soln.retcode != :Unstable || error("Solver was not stable")
  @info("Solving")
  prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
  soln = solve(prob, Tsit5())
  @info("Done")
end

begin
    frames = 300
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    msh = CairoMakie.mesh!(ax, s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
    Colorbar(fig[1,2], msh)
    CairoMakie.record(fig, "Brusselator_GPU.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
        msh.color = soln(t).U
    end
end