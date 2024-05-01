using Catlab
using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using CUDA
using CUDA.CUSPARSE
using LinearAlgebra
using SparseArrays
using Statistics
using Distributions
using Test
Point2D = Point2{Float64}
Point3D = Point3{Float64}
CUDA.allowscalar(false)

function RMSE(a, b)
  return sqrt(mean((a .- b).^2))
end

@testset "Heat Equation Float64" begin
  Heat = @decapode begin
    U::Form0
    ∂ₜ(U) == 100 * Δ(U)
  end

  tₑ = 8.5
  constants_and_parameters = ()

  s = loadmesh(Icosphere(4))
  sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(sd[:point]) do (x,_)
    return x
  end

  # CPU Setup and Solve
  cpu_sim = eval(gensim(Heat))
  cpu_fₘ = cpu_sim(sd, nothing, DiagonalHodge())

  cpu_u₀ = ComponentArray(U=U)

  prob = ODEProblem(cpu_fₘ, cpu_u₀, (0, tₑ), constants_and_parameters)
  cpu_soln = solve(prob, Tsit5(), save_everystep=false)

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Heat, code_target=cuda()))
  cuda_fₘ = cuda_sim(sd, nothing, DiagonalHodge())
  
  cuda_u₀ = ComponentArray(U=CuArray{Float64}(U))
  prob = ODEProblem(cuda_fₘ, cuda_u₀, (0, tₑ), constants_and_parameters)
  cuda_soln = solve(prob, Tsit5(), save_everystep=false)

  @test all(isapprox(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U); atol=1e-12))
  @test RMSE(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U)) < 1e-13 
end

@testset "Heat Equation Float32" begin
  Heat = @decapode begin
    U::Form0
    ∂ₜ(U) == 100 * Δ(U)
  end

  tₑ = 11.5f0
  constants_and_parameters = ()

  s = loadmesh(Icosphere(3))
  sd = EmbeddedDeltaDualComplex2D{Bool, Float32, Point3{Float32}}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(sd[:point]) do (x,_)
    return x
  end

  # CPU Setup and Solve
  cpu_sim = eval(gensim(Heat, stateeltype=Float32))
  cpu_fₘ = cpu_sim(sd, nothing, DiagonalHodge())

  cpu_u₀ = ComponentArray(U=U)

  prob = ODEProblem(cpu_fₘ, cpu_u₀, (0, tₑ), constants_and_parameters)
  cpu_soln = solve(prob, Tsit5(), save_everystep=false)

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Heat, code_target=cuda(), stateeltype=Float32))
  cuda_fₘ = cuda_sim(sd, nothing, DiagonalHodge())
  
  cuda_u₀ = ComponentArray(U=CuArray{Float32}(U))
  prob = ODEProblem(cuda_fₘ, cuda_u₀, (0, tₑ), constants_and_parameters)
  cuda_soln = solve(prob, Tsit5(), save_everystep=false)

  @test all(isapprox(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U); atol=1e-5))
  @test RMSE(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U)) < 1e-6
end

@testset "Brusselator Float64" begin
  Brusselator = @decapode begin
    (U, V)::Form0{X}
    (U2V)::Form0{X}
    (U̇, V̇)::Form0{X}

    (α)::Constant{X}
    F::Parameter{X}

    U2V == (U .* U) .* V
    aDU == (α * Δ(U))

    U̇ == 1 + U2V - (4.4 * U) + aDU + F
    V̇ == (3.4 * U) - U2V + aDU

    ∂ₜ(U) == U̇
    ∂ₜ(V) == V̇
  end

  tₑ = 25

  s = loadmesh(Icosphere(5))
  sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(sd[:point]) do (_,y,_)
    abs(y)
  end
  
  V = map(sd[:point]) do (x,_,_)
    abs(x)
  end
  
  # TODO: Try making this sparse.
  F₁ = map(sd[:point]) do (_,_,z)
    z ≥ 0.8 ? 5.0 : 0.0
  end
    
  F₂ = zeros(Float64, nv(sd))

  α = 0.001
  
  # CPU Setup and Solve
  cpu_sim = eval(gensim(Brusselator))
  cpu_fₘ = cpu_sim(sd, nothing, DiagonalHodge())

  cpu_u₀ = ComponentArray(U=U, V=V)

  cpu_constants_and_parameters = (
    α = α,
    F = t -> t ≥ 1.1 ? F₂ : F₁)  

  prob = ODEProblem(cpu_fₘ, cpu_u₀, (0, tₑ), cpu_constants_and_parameters)
  cpu_soln = solve(prob, Tsit5())

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Brusselator, code_target=cuda()))
  cuda_fₘ = cuda_sim(sd, nothing, DiagonalHodge())
  
  cuda_u₀ = ComponentArray(U=CuArray(U), V=CuArray(V))

  cuda_F₁ = CuArray(F₁)
   
  cuda_F₂ = CuArray(F₂)

  cuda_constants_and_parameters = (
    α = α,
    F = t -> t ≥ 1.1 ? cuda_F₂ : cuda_F₁)  

  prob = ODEProblem(cuda_fₘ, cuda_u₀, (0, tₑ), cuda_constants_and_parameters)
  cuda_soln = solve(prob, Tsit5())
  
  @test all(isapprox(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U); atol=1e-11))
  @test RMSE(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U)) < 1e-13
end

@testset "Klausmeier Float64" begin
  Klausmeier = @decapode begin
    (n,w)::DualForm0
    dX::Form1
    (a,ν,m)::Constant
  
    ∂ₜ(w) == a - w - w * n^2 + ν * L(dX, w)
    ∂ₜ(n) == w * n^2 - m*n + Δ(n)
  end
  
  Klausmeier[9, :type] = :DualForm0
  Klausmeier[10, :type] = :DualForm0
  Klausmeier[15, :type] = :DualForm0

  constants_and_parameters = (m = 0.45,
         a = 0.94,
         ν = 182.5)

  tₑ = 10.0

  function circle(n, c)
    s = EmbeddedDeltaSet1D{Bool, Point2D}()
    map(range(0, 2pi - (pi/(2^(n-1))); step=pi/(2^(n-1)))) do t
      add_vertex!(s, point=Point2D(cos(t),sin(t))*(c/2pi))
    end
    add_edges!(s, 1:(nv(s)-1), 2:nv(s))
    add_edge!(s, nv(s), 1)
    sd = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s)
    subdivide_duals!(sd, Circumcenter())
    s,sd
  end
  s,sd = circle(12, 500)

  n_dist = Distributions.Normal(pi)
  n = [pdf(n_dist, t)*(√(2pi))*7.2 + 0.08 - 5e-2 for t in range(0,2pi; length=ne(sd))]

  w_dist = Distributions.Normal(pi, 20)
  w = [pdf(w_dist, t) for t in range(0,2pi; length=ne(sd))]

  dX = sd[:length]

  lap_mat = dec_hodge_star(1,sd) * dec_differential(0,sd) * dec_inv_hodge_star(0,sd) * dec_dual_derivative(0,sd)

  # CPU Setup and Solve
  cpu_sim = eval(gensim(Klausmeier, dimension=1))

  function cpu_generate(sd, my_symbol; hodge=DiagonalHodge())
    op = @match my_symbol begin
      :Δ => x -> begin
        lap_mat * x
      end
    end
    return op
  end

  cpu_fₘ = cpu_sim(sd, cpu_generate, DiagonalHodge())

  cpu_u₀ = ComponentArray(n = n, w = w, dX = dX)

  prob = ODEProblem(cpu_fₘ, cpu_u₀, (0.0, tₑ), constants_and_parameters)
  cpu_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:n, :w]);

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Klausmeier, dimension=1, code_target=cuda()))
  cuda_lap_mat = CuSparseMatrixCSC(lap_mat)
  function cuda_generate(sd, my_symbol; hodge=DiagonalHodge())
    op = @match my_symbol begin
      :Δ => x -> begin
      cuda_lap_mat * x
      end
    end
    return op
  end

  cuda_fₘ = cuda_sim(sd, cuda_generate, DiagonalHodge())

  cuda_u₀ = ComponentArray(n = CuArray(n), w = CuArray(w), dX = CuArray(dX))

  prob = ODEProblem(cuda_fₘ, cuda_u₀, (0.0, tₑ), constants_and_parameters)
  cuda_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:n, :w]);

  @test all(isapprox(cpu_soln(tₑ).n, Array(cuda_soln(tₑ).n); atol=1e-11))
  @test RMSE(cpu_soln(tₑ).n, Array(cuda_soln(tₑ).n)) < 1e-13
end

@testset "Klausmeier Float32" begin
  Klausmeier = @decapode begin
    (n,w)::DualForm0
    dX::Form1
    (a,ν,m)::Constant
  
    ∂ₜ(w) == a - w - w * n^2 + ν * L(dX, w)
    ∂ₜ(n) == w * n^2 - m*n + Δ(n)
  end
  
  Klausmeier[9, :type] = :DualForm0
  Klausmeier[10, :type] = :DualForm0
  Klausmeier[15, :type] = :DualForm0

  constants_and_parameters = (m = 0.45f0,
         a = 0.94f0,
         ν = 182.5f0)

  tₑ = 3.0f0

  function circle(n, c)
    s = EmbeddedDeltaSet1D{Bool, Point2D}()
    map(range(0, 2pi - (pi/(2^(n-1))); step=pi/(2^(n-1)))) do t
      add_vertex!(s, point=Point2D(cos(t),sin(t))*(c/2pi))
    end
    add_edges!(s, 1:(nv(s)-1), 2:nv(s))
    add_edge!(s, nv(s), 1)
    sd = EmbeddedDeltaDualComplex1D{Bool, Float32, Point2{Float32}}(s)
    subdivide_duals!(sd, Circumcenter())
    s,sd
  end
  s,sd = circle(13, 500)

  n_dist = Distributions.Normal(pi)
  n = Float32.([pdf(n_dist, t)*(√(2pi))*7.2 + 0.08 - 5e-2 for t in range(0,2pi; length=ne(sd))])

  w_dist = Distributions.Normal(pi, 20)
  w = Float32.([pdf(w_dist, t) for t in range(0,2pi; length=ne(sd))])

  dX = sd[:length]

  lap_mat = dec_hodge_star(1,sd) * dec_differential(0,sd) * dec_inv_hodge_star(0,sd) * dec_dual_derivative(0,sd)

  # CPU Setup and Solve
  cpu_sim = eval(gensim(Klausmeier, dimension=1, stateeltype=Float32))

  function cpu_generate(sd, my_symbol; hodge=DiagonalHodge())
    op = @match my_symbol begin
      :Δ => x -> begin
        lap_mat * x
      end
    end
    return op
  end

  cpu_fₘ = cpu_sim(sd, cpu_generate, DiagonalHodge())

  cpu_u₀ = ComponentArray(n = n, w = w, dX = dX)

  prob = ODEProblem(cpu_fₘ, cpu_u₀, (0.0, tₑ), constants_and_parameters)
  cpu_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:n, :w]);

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Klausmeier, dimension=1, code_target=cuda(), stateeltype=Float32))
  cuda_lap_mat = CuSparseMatrixCSC(lap_mat)
  function cuda_generate(sd, my_symbol; hodge=DiagonalHodge())
    op = @match my_symbol begin
      :Δ => x -> begin
      cuda_lap_mat * x
      end
    end
    return op
  end

  cuda_fₘ = cuda_sim(sd, cuda_generate, DiagonalHodge())

  cuda_u₀ = ComponentArray(n = CuArray(n), w = CuArray(w), dX = CuArray(dX))

  prob = ODEProblem(cuda_fₘ, cuda_u₀, (0.0, tₑ), constants_and_parameters)
  cuda_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:n, :w]);

  @test RMSE(cpu_soln(tₑ).n, Array(cuda_soln(tₑ).n)) < 1e-5
end

