using ACSets
using CUDA, CUDA.CUSPARSE
using CombinatorialSpaces
using ComponentArrays
using Decapodes
using DiagrammaticEquations
using DiffEqGPU
using Distributions
using LinearAlgebra
using MLStyle
using OrdinaryDiffEqTsit5
using Statistics
using Test
CUDA.allowscalar(false)

function RMSE(a, b)
  return sqrt(mean((a .- b).^2))
end

function circle(n, c, float_type)
  s = EmbeddedDeltaSet1D{Bool, Point2d}()
  map(range(0, 2pi - (pi/(2^(n-1))); step=pi/(2^(n-1)))) do t
    add_vertex!(s, point=Point2d(cos(t),sin(t))*(c/2pi))
  end
  add_edges!(s, 1:(nv(s)-1), 2:nv(s))
  add_edge!(s, nv(s), 1)
  sd = EmbeddedDeltaDualComplex1D{Bool, float_type, Point2d}(s)
  subdivide_duals!(sd, Circumcenter())
  s, sd
end

Klausmeier = @decapode begin
  (n,w)::DualForm0
  dX::Form1
  (a,ν,m)::Constant

  ∂ₜ(w) == a - w - w * n^2 + ν * L(dX, w)
  ∂ₜ(n) == w * n^2 - m*n + Δ(n)
end

@testset "Heat Equation Float64" begin
  Heat = @decapode begin
    U::Form0
    ∂ₜ(U) == 100 * Δ(U)
  end

  tₑ = 1.0

  s = loadmesh(Icosphere(4))
  sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(p -> p[1], point(sd))

  # CPU Setup and Solve
  cpu_sim  = eval(gensim(Heat))
  cpu_fₘ   = cpu_sim(sd, nothing, DiagonalHodge())
  cpu_u₀   = ComponentArray(U=U)
  cpu_prob = ODEProblem(cpu_fₘ, cpu_u₀, (0, tₑ), ())
  cpu_soln = solve(cpu_prob, Tsit5(), save_everystep=false)

  # CUDA Setup and Solve
  cuda_sim  = eval(gensim(Heat, code_target=CUDATarget()))
  cuda_fₘ   = cuda_sim(sd, nothing, DiagonalHodge())
  cuda_u₀   = ComponentArray(U=CuArray{Float64}(U))
  cuda_prob = ODEProblem(cuda_fₘ, cuda_u₀, (0, tₑ), ())
  cuda_soln = solve(cuda_prob, Tsit5(), save_everystep=false)

  @test all(isapprox(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U); atol=1e-12))
  @test RMSE(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U)) < 1e-13
end

@testset "Heat Equation Float32" begin
  # Explicitly define the Heat equation with a Float32.
  Heat = @decapode begin
    U::Form0
    ∂ₜ(U) == 100f0 * Δ(U)
  end

  tₑ = 1.0f0

  s = loadmesh(Icosphere(3))
  sd = EmbeddedDeltaDualComplex2D{Bool, Float32, Point3{Float32}}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(p -> p[1], point(sd))

  # CPU Setup and Solve
  cpu_sim = eval(gensim(Heat, stateeltype=Float32))
  cpu_fₘ = cpu_sim(sd, nothing, DiagonalHodge())

  cpu_u₀ = ComponentArray(U=U)

  prob = ODEProblem(cpu_fₘ, cpu_u₀, (0, tₑ), ())
  cpu_soln = solve(prob, Tsit5(), save_everystep=false)

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Heat, code_target=CUDATarget(), stateeltype=Float32))
  cuda_fₘ = cuda_sim(sd, nothing, DiagonalHodge())

  cuda_u₀ = ComponentArray(U=CuArray{Float32}(U))
  prob = ODEProblem(cuda_fₘ, cuda_u₀, (0, tₑ), ())
  cuda_soln = solve(prob, Tsit5(), save_everystep=false)

  @test all(isapprox(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U); atol=1e-5))
  @test RMSE(cpu_soln(tₑ).U, Array(cuda_soln(tₑ).U)) < 1e-6
end

@testset "Brusselator Float64" begin
  Brusselator = @decapode begin
    (U, V)::Form0
    α::Constant
    F::Parameter

    U2V == (U .* U) .* V
    aDU == (α * Δ(U))

    ∂ₜ(U) == 1 + U2V - (4.4 * U) + aDU + F
    ∂ₜ(V) == (3.4 * U) - U2V + aDU
  end

  tₑ = 25

  s = loadmesh(Icosphere(5))
  sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
  subdivide_duals!(sd, Circumcenter())

  U = map(p -> abs(p[2]), point(sd))
  V = map(p -> abs(p[1]), point(sd))

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
  cpu_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:U])

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Brusselator, code_target=CUDATarget()))
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
  constants_and_parameters = (m = 0.45, a = 0.94, ν = 182.5)
  s, sd = circle(12, 500, Float64)

  n_dist = Distributions.Normal(pi)
  n = [pdf(n_dist, t)*(√(2pi))*7.2 + 0.08 - 5e-2 for t in range(0,2pi; length=ne(sd))]

  w_dist = Distributions.Normal(pi, 20)
  w = [pdf(w_dist, t) for t in range(0,2pi; length=ne(sd))]

  dX = sd[:length]

  lap_mat = dec_hodge_star(1, sd)     *
            dec_differential(0, sd)   *
            dec_inv_hodge_star(0, sd) *
            dec_dual_derivative(0, sd)

  tₑ = 2.5

  # CPU Setup and Solve
  cpu_sim = eval(gensim(Klausmeier, dimension=1))

  function cpu_generate(sd, my_symbol; hodge=DiagonalHodge())
    @match my_symbol begin
      :Δᵈ₀ => x -> lap_mat * x
      _ => error("Unmatched operator $my_symbol")
    end
  end

  cpu_fₘ   = cpu_sim(sd, cpu_generate, DiagonalHodge())
  cpu_u₀   = ComponentArray(n = n, w = w, dX = dX)
  cpu_prob = ODEProblem(cpu_fₘ, cpu_u₀, (0.0, tₑ), constants_and_parameters)
  cpu_soln = solve(cpu_prob, Tsit5(), save_everystep=false, save_idxs=[:n]);

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Klausmeier, dimension=1, code_target=CUDATarget()))
  cuda_lap_mat = CuSparseMatrixCSC(lap_mat)
  function cuda_generate(sd, my_symbol; hodge=DiagonalHodge())
    @match my_symbol begin
      :Δᵈ₀ => x -> cuda_lap_mat * x
      _ => error("Unmatched operator $my_symbol")
    end
  end

  cuda_fₘ   = cuda_sim(sd, cuda_generate, DiagonalHodge())
  cuda_u₀   = ComponentArray(n = CuArray(n), w = CuArray(w), dX = CuArray(dX))
  cuda_prob = ODEProblem(cuda_fₘ, cuda_u₀, (0.0, tₑ), constants_and_parameters)
  cuda_soln = solve(cuda_prob, Tsit5(), save_everystep=false, save_idxs=[:n]);

  @test all(isapprox(cpu_soln(tₑ).n, Array(cuda_soln(tₑ).n); atol=1e-11))
  @test RMSE(cpu_soln(tₑ).n, Array(cuda_soln(tₑ).n)) < 1e-13
end

@testset "Klausmeier Float32" begin
  constants_and_parameters = (m = 0.45f0, a = 0.94f0, ν = 182.5f0)
  s, sd = circle(13, 500, Float32)

  n_dist = Distributions.Normal(pi)
  n = map(t -> Float32(pdf(n_dist, t)*(√(2pi))*7.2 + 0.08 - 5e-2),
    range(0, 2pi; length=ne(sd)))

  w_dist = Distributions.Normal(pi, 20)
  w = map(t -> Float32(pdf(w_dist, t)), range(0, 2pi; length=ne(sd)))

  dX = sd[:length]

  lap_mat = dec_hodge_star(1, sd)     *
            dec_differential(0, sd)   *
            dec_inv_hodge_star(0, sd) *
            dec_dual_derivative(0, sd)

  tₑ = 2.5f0

  # CPU Setup and Solve
  cpu_sim  = eval(gensim(Klausmeier, dimension=1, stateeltype=Float32))
  cpu_fₘ   = cpu_sim(sd, nothing, DiagonalHodge())
  cpu_u₀   = ComponentArray(n = n, w = w, dX = dX)
  prob     = ODEProblem(cpu_fₘ, cpu_u₀, (0.0, tₑ), constants_and_parameters)
  cpu_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:n]);

  # CUDA Setup and Solve
  cuda_sim = eval(gensim(Klausmeier, dimension=1, code_target=CUDATarget(), stateeltype=Float32))
  cuda_lap_mat = CuSparseMatrixCSC(lap_mat)
  function cuda_generate(sd, my_symbol; hodge=DiagonalHodge())
    @match my_symbol begin
      :Δᵈ₀ => x -> cuda_lap_mat * x
      _ => error("Unmatched operator $my_symbol")
    end
  end

  cuda_fₘ   = cuda_sim(sd, cuda_generate, DiagonalHodge())
  cuda_u₀   = ComponentArray(n = CuArray(n), w = CuArray(w), dX = CuArray(dX))
  prob      = ODEProblem(cuda_fₘ, cuda_u₀, (0.0, tₑ), constants_and_parameters)
  cuda_soln = solve(prob, Tsit5(), save_everystep=false, save_idxs=[:n]);

  @test RMSE(cpu_soln(tₑ).n, Array(cuda_soln(tₑ).n)) < 1e-5
end

@testset "Ensemble GPU Sims" begin
  # Define model.
  Heat = @decapode begin
    C::Form0
    D::Constant
    ∂ₜ(C) == D*Δ(C)
  end
  # Define domain.
  s, sd = circle(7, 500, Float64)
  # Create initial data.
  Csin = map(p -> sin(p[1]), point(s))
  Ccos = map(p -> cos(p[1]), point(s))
  constants_and_parameters = (D = 0.001,)

  # In this test, set up a Decapodes GPU simulation and parallelize with usual threads.
  sim = eval(gensim(Heat, dimension=1, code_target=CUDATarget()))
  fₘ = sim(sd, nothing)
  tₑ = 1e8
  C = CuArray(stack([Csin, Ccos]))
  u₀ = ComponentArray(C=Csin,)
  ode_prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
  ens_prob = EnsembleProblem(ode_prob,
    prob_func = (prob, i, repeat) ->
      remake(prob, u0=ComponentArray(C=C[:,i])))
  soln = solve(ens_prob, Tsit5(), EnsembleThreads(); trajectories=2)
  @test all(Array(soln[1].u[1]) .== Csin)
  @test all(Array(soln[1].u[1]) .!= Ccos)
  @test all(Array(soln[2].u[1]) .!= Csin)
  @test all(Array(soln[2].u[1]) .== Ccos)

  # In this test, set up a Decapodes CPU simulation and parallelize an ensemble GPU array.
  # This fails, because the SciML GPU compiler cannot compile the sparse matrix
  # operation, as it is not isbitstype.
  #=
  ERROR: GPU compilation of MethodInstance for DiffEqGPU.gpu_gpu_kernel(::KernelAbstractions.CompilerMetadata{…}, ::var"#f#314"{…}, ::CuDeviceMatrix{…}, ::CuDeviceMatrix{…}, ::CuDeviceMatrix{…}, ::Float64) failed
  KernelError: passing non-bitstype argument
  
  Argument 3 to your kernel function is of type var"#f#314"{var"#311#313"{SparseArrays.SparseMatrixCSC{Float64, Int32}}}, which is not a bitstype:
    .GenSim-ConMat_0 is of type var"#311#313"{SparseArrays.SparseMatrixCSC{Float64, Int32}} which is not isbits.
      .GenSim-M_GenSim-ConMat_0 is of type SparseArrays.SparseMatrixCSC{Float64, Int32} which is not isbits.
  =#
  @test_broken let
    sim = eval(gensim(Heat, dimension=1, preallocate=false))
    fₘ = sim(sd, nothing)
    tₑ = 1e8
    C = stack([Csin, Ccos])
    u₀ = ComponentArray(C=Csin,)
    ode_prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
    ens_prob = EnsembleProblem(ode_prob,
      prob_func = (prob, i, repeat) ->
        remake(prob, u0=ComponentArray(C=C[:,i])))
    soln = solve(ens_prob, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()); trajectories=2)
    @test all(Array(soln[1].u[1]) .== Csin)
    @test all(Array(soln[1].u[1]) .!= Ccos)
    @test all(Array(soln[2].u[1]) .!= Csin)
    @test all(Array(soln[2].u[1]) .== Ccos)
    true
  end
end
