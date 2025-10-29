using ACSets
using Catlab
using CombinatorialSpaces
using ComponentArrays
using Decapodes
using DiagrammaticEquations
using Distributions
using GeometryBasics: Point2, Point3
using LinearAlgebra
using MLStyle
using OrdinaryDiffEqTsit5
using Test
using Random
Point3D = Point3{Float64}

import Decapodes: default_dec_matrix_generate

flatten(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function test_hodge(k, sd::HasDeltaSet, hodge)
  hodge = ⋆(k,sd,hodge=hodge)
  x-> hodge * x
end

function test_inverse_hodge(k, sd::HasDeltaSet, hodge)
  invhodge = inv_hodge_star(k,sd,hodge)
  x-> invhodge * x
end

function test_differential(k, sd::HasDeltaSet)
  diff = d(k,sd)
  x-> diff * x
end

function test_dual_differential(k, sd::HasDeltaSet)
  dualdiff = dual_derivative(k,sd)
  x-> dualdiff * x
end

function test_codifferential(k, sd::HasDeltaSet, hodge)
  codiff = δ(k, sd, hodge, nothing)
  x -> codiff * x
end

function test_laplace_de_rham(k, sd::HasDeltaSet)
  lpdr = Δ(k, sd)
  x -> lpdr * x
end

function dec_laplace_beltrami(k, sd::HasDeltaSet)
  lpbt = ∇²(k, sd)
  x -> lpbt * x
end

function generate(sd, my_symbol)
  op = @match my_symbol begin
    :⋆₁ => test_hodge(1, sd, DiagonalHodge())
    :⋆₀⁻¹ => test_inverse_hodge(0, sd, DiagonalHodge())
    :dual_d₁ => test_dual_differential(1, sd)
    _ => default_dec_generate_2D(sd, my_symbol)
  end
  return (args...) ->  op(args...)
end

@testset "Simulation Generation" begin
# Mesh and ICs to use for these tests:
torus = loadmesh(Torus_30x10())
c_dist = MvNormal([5, 5], LinearAlgebra.Diagonal(map(abs2, [1.5, 1.5])))
c = [pdf(c_dist, [p[1], p[2]]) for p in torus[:point]]
u₀ = ComponentArray(C=c)
du = ComponentArray(C=zero(c))

# Three Decapodes variations, with k as constant, parameter, or literal.
DiffusionWithConstant = @decapode begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}
  k::Constant{ℝ}

  # Fick's first law
  ϕ == k * d₀(C)
  # Diffusion equation
  Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
  ∂ₜ(C) == Ċ
end

DiffusionWithParameter = @decapode begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}
  k::Parameter{ℝ}

  # Fick's first law
  ϕ == k * d₀(C)
  # Diffusion equation
  Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
  ∂ₜ(C) == Ċ
end

DiffusionWithLiteral = @decapode begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ == 3 * d₀(C)
  # Diffusion equation
  Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
  ∂ₜ(C) == Ċ
end

# Verify the variable accessors.
@test Decapodes.get_vars_code(DiffusionWithConstant, [:k], Float64, CPUTarget()).args[2] == :(k = __p__.k)
@test infer_state_names(DiffusionWithConstant) == [:C, :k]

@test infer_state_names(DiffusionWithParameter) == [:C, :k]
@test Decapodes.get_vars_code(DiffusionWithParameter, [:k], Float64, CPUTarget()).args[2] == :(k = __p__.k(__t__))

@test infer_state_names(DiffusionWithLiteral) == [:C]
# TODO: Fix proper Expr equality, the Float64 does not equate here
# @test Decapodes.get_vars_code(DiffusionWithLiteral, [Symbol("3")]).args[2] == :(var"3"::Float64 = 3.0)
@test Decapodes.get_vars_code(DiffusionWithLiteral, [Symbol("3")], Float64, CPUTarget()).args[2].args[1] == :(var"3")
@test Decapodes.get_vars_code(DiffusionWithLiteral, [Symbol("3")], Float64, CPUTarget()).args[2].args[2] == 3.0

# Test that simulations generated from these return the same result.
f = evalsim(DiffusionWithConstant)
f_with_constant = f(torus, generate)

f = evalsim(DiffusionWithParameter)
f_with_parameter = f(torus, generate)

f = eval(gensim(expand_operators(DiffusionWithLiteral)))
f_with_literal = f(torus, generate)

f_with_constant(du, u₀, (k=3.0,), 0)
fc_res = copy(du.C)
f_with_parameter(du, u₀, (k=t->3.0,), 0)
fp_res = copy(du.C)
f_with_literal(du, u₀, NamedTuple(), 0)
fl_res = copy(du.C)

@test norm(fc_res - fp_res) < 1e-4
@test norm(fc_res - fl_res) < 1e-4

# Test same but with no preallocating

f = evalsim(DiffusionWithLiteral, preallocate=false)
f_noalloc = f(torus, generate)

f = evalsim(DiffusionWithLiteral)
f_alloc = f(torus, generate)

f_noalloc(du, u₀, NamedTuple(), 0)
f_nal = copy(du.C)
f_alloc(du, u₀, NamedTuple(), 0)
f_al = copy(du.C)

@test f_nal == f_al

end

# Testing done based on the original gensim
# -----------------------------------------

# Testing Brusselator
@testset "Brusselator Simulation" begin

  function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :Δ₀ => test_laplace_de_rham(0, sd)
      _ => default_dec_generate_2D(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
  end

  begin
    primal_earth = loadmesh(Icosphere(1))
    nploc = argmax(x -> x[3], primal_earth[:point])
    primal_earth[:edge_orientation] = false
    orient!(primal_earth)
    earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
    subdivide_duals!(earth, Circumcenter());
  end

  Brusselator = SummationDecapode(parse_decapode(quote
      (U, V)::Form0{X}
      (U2V, aTU)::Form0{X}
      (U̇, V̇)::Form0{X}

      (α, One)::Constant{X}
      (F)::Parameter{X}

      U2V == (U .* U) .* V
      aTU == α * Δ(U)

      U̇ == One + U2V - (4.4 * U) + aTU + F
      V̇ == (3.4 * U) - U2V + aTU

      ∂ₜ(U) == U̇
      ∂ₜ(V) == V̇
  end))

  function old_simulate(mesh, operators)
    begin
        Δ₀ = generate(mesh, :Δ₀)
        # (.*) = operators(mesh, :.*)
    end
    return begin
            f(du, u, p, t) = begin
                    begin
                        U = u.U
                        V = u.V
                        α = p.α
                        One = p.One
                        F = p.F(t)
                        var"4.4" = 4.4
                        var"3.4" = 3.4
                    end
                    var"•4" = Δ₀(U)
                    var"•2" = U .* U
                    U2V = var"•2" .* V
                    aTU = α * var"•4"
                    var"•6" = var"4.4" * U
                    var"•3" = var"3.4" * U
                    var"•1" = var"•3" - U2V
                    sum_1 = One + U2V
                    V̇ = var"•1" + aTU
                    var"•5" = sum_1 - var"•6"
                    U̇ = var"•5" + aTU + F
                    du .= 0.0
                    begin
                        (du.U) .= U̇
                        (du.V) .= V̇
                    end
                end
        end
  end

  fₙ = old_simulate(earth, generate)

  new_sim = evalsim(Brusselator)
  fₘ = new_sim(earth, generate)

  begin
    U = map(earth[:point]) do (_,y,_)
      abs(y)
    end

    V = map(earth[:point]) do (x,_,_)
      abs(x)
    end

    One = ones(nv(earth))

    F₁ = map(earth[:point]) do (_,_,z)
      z ≥ 0.8 ? 5.0 : 0.0
    end

    F₂ = zeros(nv(earth))

    constants_and_parameters = (
        α = 0.001,
        F = t -> t ≥ 1.1 ? F₁ : F₂,
        One = One)

    u₀ = ComponentArray(U=U,V=V)
    tₑ = 11.5
    prob = ODEProblem(fₙ,u₀,(0, tₑ), constants_and_parameters)
    old_soln = solve(prob, Tsit5())
  end

  begin
      U = map(earth[:point]) do (_,y,_)
        abs(y)
      end

      V = map(earth[:point]) do (x,_,_)
        abs(x)
      end

      One = ones(nv(earth))

      F₁ = map(earth[:point]) do (_,_,z)
        z ≥ 0.8 ? 5.0 : 0.0
      end

      F₂ = zeros(nv(earth))

      constants_and_parameters = (
          α = 0.001,
          F = t -> t ≥ 1.1 ? F₁ : F₂,
          One = One)

      u₀ = ComponentArray(U=U,V=V)
      tₑ = 11.5
      prob = ODEProblem(fₘ,u₀,(0, tₑ), constants_and_parameters)
      new_soln = solve(prob, Tsit5())
  end

  @test old_soln.u ≈ new_soln.u
end

# Testing Budyko-Sellers
@testset "Budyko-Sellers Simulation" begin
  # This is a 1D DEC test.
  # The dimension impacts the allocation of DualForms.
  budyko_sellers = @decapode begin
    (Q,Tₛ)::Form0
    (α,A,B,C,D,cosϕᵖ,cosϕᵈ)::Constant

    Tₛ̇ == ∂ₜ(Tₛ)
    ASR == (1 .- α) .* Q
    OLR == A .+ (B .* Tₛ)
    HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))

    Tₛ̇ == (ASR - OLR + HT) ./ C
  end
  infer_types!(budyko_sellers, dim=1)
  resolve_overloads!(budyko_sellers, dim=1)

  # This test ensures that the next one does not break, since it depends on
  # arbitrary internal variable naming.
  @test budyko_sellers[only(incident(budyko_sellers, Symbol("•1"), :name)), :type] == :DualForm0
  # A dual 0-form consists of ne(s) floats.
  @test occursin("var\"__•1\" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))",
    repr(gensim(budyko_sellers, dimension=1)))
end

@testset "Gensim Transformations" begin

  function count_contractions(e::Expr)
    block = e.args[2].args[2].args[5]
    length(block.args) - 1
  end

  count_contractions(d::SummationDecapode) = count_contractions(gensim(d))

  begin
    primal_earth = loadmesh(Icosphere(1))
    nploc = argmax(x -> x[3], primal_earth[:point])
    primal_earth[:edge_orientation] = false
    orient!(primal_earth)
    earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
    subdivide_duals!(earth, Circumcenter());
  end

  @test_throws "Unmatched operator Test" default_dec_generate(earth, Symbol("Test"))

  # Testing transforming negation into multiplication by -1
  neg_transform = @decapode begin
    (A,B)::Form0

    B == ∂ₜ(A)

    B == -(neg(-1.0 * A))
  end

  sim = eval(gensim(neg_transform))
  f = sim(earth, default_dec_generate)
  A = ones(nv(earth))
  u = ComponentArray(A=A)
  constants_and_parameters = ()
  f(u, u, constants_and_parameters, 0)

  @test u.A == -1 * ones(nv(earth))


  # Testing simple contract operations
  single_contract = @decapode begin
    (A,C,E)::Form0
    (D,F)::Form2

    B == ∂ₜ(A)
    D == ∂ₜ(C)

    B == ⋆(⋆(A))
    D == d(d(C))
    F == d(d(E))
  end
  @test 4 == count_contractions(single_contract)

  @test 0 == count_contractions(gensim(single_contract; contract=false))

  f = gensim(single_contract)
  @test f.args[2].args[2].args[5].args[[2,4]] == [
    :(var"GenSim-M_GenSim-ConMat_0" = var"GenSim-M_d₁" * var"GenSim-M_d₀"),
    :(var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_⋆₀")]

  sim = eval(gensim(single_contract))
  f = sim(earth, default_dec_generate)
  A = 2 * ones(nv(earth))
  C = ones(nv(earth))
  E = ones(nv(earth))
  u = ComponentArray(A=A, C=C, E=E)
  du = ComponentArray(A=zeros(nv(earth)), C=zeros(ntriangles(earth)), E=zeros(ntriangles(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A ≈ 2 * ones(nv(earth))
  @test du.C == zeros(ntriangles(earth))
  @test du.E == zeros(ntriangles(earth))

  # Testing contraction interrupted by summation
  contract_with_summation = @decapode begin
    (A)::Form0
    (D)::Form2

    C == ∂ₜ(E)
    D == ∂ₜ(A)

    B == ⋆(⋆(A))
    C == B + B

    D == d(d(C))
  end
  @test 4 == count_contractions(contract_with_summation)

  sim = eval(gensim(contract_with_summation))
  f = sim(earth, default_dec_generate)
  A = 2 * ones(nv(earth))
  E_dec = ones(nv(earth))
  u = ComponentArray(A=A, E=E_dec)
  du = ComponentArray(A=zeros(ntriangles(earth)), E=zeros(nv(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == zeros(ntriangles(earth))
  @test du.E ≈ 4 * ones(nv(earth))

  # Testing contraction interrupted by op2
  contract_with_op2 = @decapode begin
    (A)::Form0
    (D)::Form2

    C == ∂ₜ(E)
    D == ∂ₜ(A)

    B == ⋆(⋆(A))
    C == B * B

    D == d(d(C))
  end
  @test 4 == count_contractions(contract_with_op2)

  sim = eval(gensim(contract_with_op2, preallocate = false))
  f = sim(earth, default_dec_generate)
  A = 3 * ones(nv(earth))
  E_dec = ones(nv(earth))
  u = ComponentArray(A=A, E=E_dec)
  du = ComponentArray(A=zeros(ntriangles(earth)), E=zeros(nv(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == zeros(ntriangles(earth))
  @test du.E ≈ 9 * ones(nv(earth))

  sim = eval(gensim(contract_with_op2, preallocate = true))
  f = sim(earth, default_dec_generate)
  A = 3 * ones(nv(earth))
  E_dec = ones(nv(earth))
  u = ComponentArray(A=A, E=E_dec)
  du = ComponentArray(A=zeros(ntriangles(earth)), E=zeros(nv(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == zeros(ntriangles(earth))
  @test du.E ≈ 9 * ones(nv(earth))


  # Testing contract lines beyond the initial value
  later_contraction = @decapode begin
    (A)::Form0

    D == ∂ₜ(A)

    B == A * A
    D == ⋆(⋆(B))
  end
  @test 2 == count_contractions(later_contraction)

  sim = eval(gensim(later_contraction))
  f = sim(earth, default_dec_generate)
  A = 4 * ones(nv(earth))
  u = ComponentArray(A=A)
  du = ComponentArray(A=zeros(nv(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A ≈ 16 * ones(nv(earth))

  # Testing no contraction of single operators
  no_contraction = @decapode begin
    (A)::Form0
    (D)::Form1

    D == ∂ₜ(A)
    D == d(A)
  end
  @test 0 == count_contractions(no_contraction)

  sim = eval(gensim(no_contraction))
  f = sim(earth, default_dec_generate)
  A = [i for i in 1:nv(earth)]
  u = ComponentArray(A=A)
  du = ComponentArray(A=zeros(ne(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == CombinatorialSpaces.DiscreteExteriorCalculus.d(0, earth) * A

  # Testing no contraction of unallowed operators
  no_unallowed = @decapode begin
    (A)::Form0
    (D)::Form1

    D == ∂ₜ(A)
    D == d(k(A))
  end
  @test 0 == count_contractions(no_unallowed)

  sim = eval(gensim(no_unallowed))

  function generate_no_unallowed(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :k => (x -> 20 * x)
    end
    op
  end

  f = sim(earth, generate_no_unallowed)
  A = [i for i in 1:nv(earth)]
  u = ComponentArray(A=A)
  du = ComponentArray(A=zeros(ne(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == CombinatorialSpaces.DiscreteExteriorCalculus.d(0, earth) * 20 * A

  # Testing wedge 01 operators function
  wedges01 = @decapode begin
    (A, B)::Form0
    (C, D, E)::Form1

    D == ∂ₜ(A)
    E == ∂ₜ(B)
    F == ∂ₜ(C)


    D == (A ∧ B) ∧ C
    E == A ∧ (B ∧ C)

    F == A ∧ (C ∧ B)
  end

  sim = eval(gensim(wedges01, preallocate=false))
  f = sim(earth, default_dec_generate)
  A = ones(nv(earth))
  B = 2 * ones(nv(earth))
  C = 3 * ones(ne(earth))
  u = ComponentArray(A=A, B=B, C=C)
  du = ComponentArray(A=zeros(ne(earth)), B=zeros(ne(earth)), C=zeros(ne(earth)))

  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == du.B == du.C

  sim = eval(gensim(wedges01, preallocate=true))
  f = sim(earth, default_dec_generate)
  A = ones(nv(earth))
  B = 2 * ones(nv(earth))
  C = 3 * ones(ne(earth))
  u = ComponentArray(A=A, B=B, C=C)
  du = ComponentArray(A=zeros(ne(earth)), B=zeros(ne(earth)), C=zeros(ne(earth)))

  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == du.B == du.C

  # Testing wedge 11 operators function
  wedges11 = @decapode begin
    (A, B)::Form1
    (D, E)::Form2

    D == ∂ₜ(A)
    E == ∂ₜ(B)

    D == A ∧ B
    E == B ∧ A
  end

  sim = eval(gensim(wedges11))

  f = sim(earth, default_dec_generate)
  A = ones(ne(earth))
  B = ones(ne(earth))
  u = ComponentArray(A=A, B=B)
  du = ComponentArray(A=zeros(ntriangles(earth)), B=zeros(ntriangles(earth)))

  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test all(isapprox.(du.A, zeros(ntriangles(earth)); atol = 1e-15))
  @test all(isapprox.(du.B, zeros(ntriangles(earth)); atol = 1e-15))
  @test all(isapprox.(du.A, du.B; atol = 1e-15))

  # Testing wedge 02 operators function
  wedges02 = @decapode begin
    A::Form0
    B::Form2
    (D, E)::Form2

    D == ∂ₜ(A)
    E == ∂ₜ(B)

    D == A ∧ B
    E == B ∧ A
  end

  sim = eval(gensim(wedges02))

  f = sim(earth, default_dec_generate)
  A = ones(nv(earth))
  B = ones(ntriangles(earth))
  u = ComponentArray(A=A, B=B)
  du = ComponentArray(A=zeros(ntriangles(earth)), B=zeros(ntriangles(earth)))

  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test all(isapprox.(du.A, ones(ntriangles(earth))))
  @test all(isapprox.(du.B, ones(ntriangles(earth))))
  @test all(isapprox.(du.A, du.B))

  # Testing Geo inverse hodge 1
  GeoInvHodge1 = @decapode begin
    A::DualForm1

    B == ∂ₜ(A)
    B == ⋆(⋆(A))
  end

  sim = eval(gensim(GeoInvHodge1))

  f = sim(earth, default_dec_generate)
  A = ones(ne(earth))
  u = ComponentArray(A=A)
  du = ComponentArray(A=zeros(ne(earth)))

  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test all(isapprox.(du.A, -1 * ones(ne(earth))))

  # Testing Diagonal inverse hodge 1
  DiagonalInvHodge1 = @decapode begin
    A::DualForm1

    B == ∂ₜ(A)
    B == ⋆(⋆(A))
  end

  sim = eval(gensim(DiagonalInvHodge1))

  f = sim(earth, default_dec_generate, DiagonalHodge())
  A = ones(ne(earth))
  u = ComponentArray(A=A)
  du = ComponentArray(A=zeros(ne(earth)))

  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test all(isapprox.(du.A, -1 * ones(ne(earth))))

end

filter_lnn(arr::AbstractVector) = filter(x -> !(x isa LineNumberNode), arr)

@testset "1-D Mat Generation" begin
  Point2D = Point2{Float64}
  function generate_dual_mesh(s::HasDeltaSet1D)
    orient!(s)
    sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point2D}(s)
    subdivide_duals!(sd, Barycenter())
    sd
  end

  primal_line = EmbeddedDeltaSet1D{Bool,Point2D}()
  add_vertices!(primal_line, 3, point=[Point2D(1,0), Point2D(0,0), Point2D(0,2)])
  add_edges!(primal_line, [1,2], [2,3])
  line = generate_dual_mesh(primal_line)

  # Testing Diagonal inverse hodge 1
  DiagonalInvHodge1 = @decapode begin
    A::DualForm1

    B == ∂ₜ(A)
    B == ⋆(⋆(A))
  end
  g = gensim(DiagonalInvHodge1)
  @test g.args[2].args[2].args[3].args[2].args[2].args[3].value == :⋆₁⁻¹
  @test length(filter_lnn(g.args[2].args[2].args[3].args)) == 2
  sim = eval(g)

  # TODO: Error is being thrown here
  @test f = sim(line, default_dec_generate, DiagonalHodge()) isa Any
end

@testset "GenSim Compilation" begin

  rect = triangulated_grid(100, 100, 50, 50, Point3{Float64})
  d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}}(rect)
  subdivide_duals!(d_rect, Circumcenter())

  function generate(sd, my_symbol, hodge)
    op = @match my_symbol begin
      end
    return op
  end

  # tests that there is no variable shadowing for du, u, p, and t
  NoShadow = @decapode begin
      u::Form0
      du::Form1
      v::Form0
      p::Constant
      q::Constant
      t::Constant
  end
  symsim = gensim(NoShadow)
  sim_NS = eval(symsim)
  @test sim_NS(d_rect, generate, DiagonalHodge()) isa Any

  HeatTransfer = @decapode begin
    (HT, Tₛ)::Form0
    (D, cosϕᵖ, cosϕᵈ)::Constant
    HT == (D ./ cosϕᵖ) .* (⋆)(d(cosϕᵈ .* (⋆)(d(Tₛ))))
  end

  sim_HT = evalsim(HeatTransfer)
  @test sim_HT(d_rect, generate, DiagonalHodge()) isa Any

  Jordan_Kinderlehrer_Otto = @decapode begin
    (ρ, Ψ)::Form0
    β⁻¹::Constant
    ∂ₜ(ρ) == (∘(⋆, d, ⋆))(d(Ψ) ∧ ρ) + β⁻¹ * Δ(ρ)
  end

  sim_JKO = evalsim(Jordan_Kinderlehrer_Otto)
  @test sim_JKO(d_rect, generate, DiagonalHodge()) isa Any

  Schoedinger = @decapode begin
    (i, h, m)::Constant
    V::Parameter
    Ψ::Form0
    ∂ₜ(Ψ) == (((-1 * h ^ 2) / (2m)) * Δ(Ψ) + V * Ψ) / (i * h)
  end
  sim_Schoedinger = evalsim(Schoedinger)
  @test sim_Schoedinger(d_rect, generate, DiagonalHodge()) isa Any

  Gray_Scott = @decapode begin
    (U, V)::Form0
    UV2::Form0
    (U̇, V̇)::Form0
    (f, k, rᵤ, rᵥ)::Constant
    UV2 == U .* (V .* V)
    U̇ == (rᵤ * Δ(U) - UV2) + f * (1 .- U)
    V̇ == (rᵥ * Δ(V) + UV2) - (f + k) .* V
    ∂ₜ(U) == U̇
    ∂ₜ(V) == V̇
  end

  sim_GS = evalsim(Gray_Scott)
  @test sim_GS(d_rect, generate, DiagonalHodge()) isa Any

  Lejeune = @decapode begin
    ρ::Form0
    (μ, Λ, L)::Constant
    ∂ₜ(ρ) == (ρ * (((1 - μ) + (Λ - 1) * ρ) - ρ * ρ) + 0.5 * (L * L - ρ) * Δ(ρ)) - 0.125 * ρ * Δ(ρ) * Δ(ρ)
  end

  sim_LJ = evalsim(Lejeune)
  @test sim_LJ(d_rect, generate, DiagonalHodge()) isa Any

  Tracer = @decapode begin
    (c, C, F, c_up)::Form0
    (v, V, q)::Form1
    c_up == (((-1 * (⋆)(L(v, (⋆)(c))) - (⋆)(L(V, (⋆)(c)))) - (⋆)(L(v, (⋆)(C)))) - (∘(⋆, d, ⋆))(q)) + F
  end

  sim_Tracer = evalsim(Tracer)
  @test sim_Tracer(d_rect, generate, DiagonalHodge()) isa Any

  # Test for Halfar
  halfar_eq2 = @decapode begin
    h::Form0
    Γ::Form1
    n::Constant

    ḣ == ∂ₜ(h)
    ḣ == ∘(⋆, d, ⋆)(Γ * d(h) ∧₁₀ ((mag(♯ᵖᵖ(d(h)))^(n-1)) ∧₀₀ h^(n+2)))
  end

  glens_law = @decapode begin
    Γ::Form1
    A::Form1
    (ρ,g,n)::Constant

    Γ == (2/(n+2))*A*(ρ*g)^n
  end

  ice_dynamics_composition_diagram = @relation () begin
    dynamics(Γ,n)
    stress(Γ,n)
  end

  ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
    [Open(halfar_eq2, [:Γ,:n]),
     Open(glens_law, [:Γ,:n])])
  halfar = apex(ice_dynamics_cospan)

  resolve_overloads!(infer_types!(halfar))

  function halfar_generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :norm => x -> norm.(x)
      x => error("Unmatched operator $my_symbol")
    end
    return op
  end

  sim_Halfar = evalsim(halfar)
  @test sim_Halfar(d_rect, halfar_generate, DiagonalHodge()) isa Any

  # Test for Poisson
  eq11_inviscid_poisson = @decapode begin
    d𝐮::DualForm2
    𝐮::DualForm1
    ψ::Form0

    ψ == Δ₀⁻¹(⋆(d𝐮))
    𝐮 == ⋆(d(ψ))

    ∂ₜ(d𝐮) ==  (-1) * ∘(♭♯, ⋆₁, d̃₁)(∧ᵈᵖ₁₀(𝐮, ⋆(d𝐮)))
  end

  function poisson_generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :♭♯ => x -> nothing
      x => error("Unmatched operator $my_symbol")
    end
    return op
  end

  sim_Poisson = evalsim(eq11_inviscid_poisson)
  @test sim_Poisson(d_rect, poisson_generate, DiagonalHodge()) isa Any

  # Test for Halmo
  eq10forN2 = @decapode begin
    (𝐮,w)::DualForm1
    (P, 𝑝ᵈ)::DualForm0
    μ::Constant

    𝑝ᵈ == P + 0.5 * ι₁₁(w,w)

    ∂ₜ(𝐮) == μ * ∘(d, ⋆, d, ⋆)(w) + (-1)*⋆₁⁻¹(∧ᵈᵖ₁₀(w, ⋆(d(w)))) + d(𝑝ᵈ)
  end

  halfar_eq2 = @decapode begin
    h::Form0
    Γ::Form1
    n::Constant

    ∂ₜ(h) == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯ᵖᵖ(d(h)))^(n-1)) * avg₀₁(h^(n+2)))
  end

  glens_law = @decapode begin
    Γ::Form1
    (A,ρ,g,n)::Constant

    Γ == (2/(n+2))*A*(ρ*g)^n
  end

  ice_dynamics_composition_diagram = @relation () begin
    dynamics(Γ,n)
    stress(Γ,n)
  end

  ice_dynamics = apex(oapply(ice_dynamics_composition_diagram,
    [Open(halfar_eq2, [:Γ,:n]),
     Open(glens_law, [:Γ,:n])]))

     ice_water_composition_diagram = @relation () begin
     glacier_dynamics(ice_thickness)
     water_dynamics(flow, flow_after)

     interaction(ice_thickness, flow, flow_after)
   end

   blocking = @decapode begin
    h::Form0
    (𝐮,w)::DualForm1

    w == (1-σ(h)) ∧ᵖᵈ₀₁ 𝐮
  end

  ice_water = apex(oapply(ice_water_composition_diagram,
  [Open(ice_dynamics, [:dynamics_h]),
   Open(eq10forN2,    [:𝐮, :w]),
   Open(blocking,     [:h, :𝐮, :w])]))

  function halmo_generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :σ => x -> nothing
      :norm => x -> nothing
      _ => error("Unmatched operator $my_symbol")
    end
    return op
  end

  resolve_overloads!(infer_types!(ice_water))

  sim_Halmo = evalsim(ice_water)
  @test sim_Halmo(d_rect, halmo_generate, DiagonalHodge()) isa Any

end

@testset "Multigrid" begin
  s = triangulated_grid(1,1,1/4,sqrt(3)/2*1/4,Point3D)

  series = PrimalGeometricMapSeries(s, binary_subdivision_map, 4);

  our_mesh = finest_mesh(series)
  lap = ∇²(0,our_mesh);

  Random.seed!(1337)
  b = lap*rand(nv(our_mesh));

  inv_lap = @decapode begin
    U::Form0
    ∂ₜ(U) == Δ₀⁻¹(U)
  end

  function generate(fs, my_symbol; hodge=DiagonalHodge())
    op = @match my_symbol begin
      _ => default_dec_matrix_generate(fs, my_symbol, hodge)
    end
  end

  sim = eval(gensim(inv_lap))
  sim_mg = eval(gensim(inv_lap; multigrid=true))

  f = sim(our_mesh, generate);
  f_mg = sim_mg(series, generate);

  u = ComponentArray(U=b)
  du = similar(u)

  # Regular mesh
  f(du, u, 0, ())
  @test norm(lap*du.U-b)/norm(b) < 1e-15

  # Multigrid
  f_mg(du, u, 0, ())
  @test norm(lap*du.U-b)/norm(b) < 1e-6
end

@testset "Allocations" begin
# Test the heat equation Decapode has expected memory allocation.

  Heat = @decapode begin
    C::Form0
    D::Constant
    ∂ₜ(C) == D*Δ(C)
  end

  s = loadmesh(Icosphere(1))
  sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
  subdivide_duals!(sd, Circumcenter())
  u₀ = ComponentArray(C = map(x -> x[3], point(sd)))
  p = (D=1e-1,)
  du = copy(u₀)

  # Allocating
  sim = eval(gensim(Heat, preallocate=false))
  f = sim(sd,nothing)
  # The first call to the function makes many allocations.
  _ = @allocations f(du, u₀, p, (0,1.0)) # 55259
  _ = @allocated f(du, u₀, p, (0,1.0)) # 3962696
  # Test that subsequent calls make a reasonable amount.
  nallocs = @allocations f(du, u₀, p, (0,1.0))
  bytes = @allocated f(du, u₀, p, (0,1.0))

  @test (nallocs, bytes) <= (7, 400)

  # Not allocating
  Heat = @decapode begin
    C::Form0
    D::Constant
    ∂ₜ(C) == D*Δ(C)
  end

  sim = eval(gensim(Heat, preallocate=true))
  f = sim(sd,nothing)
  # The first call to the function makes many allocations.
  _ = @allocations f(du, u₀, p, (0,1.0)) # 55259
  _ = @allocated f(du, u₀, p, (0,1.0)) # 3962696
  # Test that subsequent calls make a reasonable amount.
  nallocs = @allocations f(du, u₀, p, (0,1.0))
  bytes = @allocated f(du, u₀, p, (0,1.0))

  @test (nallocs, bytes) <= (6, 80)

end

@testset "Ensemble Simulations" begin
  # Define Model
  Heat = @decapode begin
    C::Form0
    D::Constant
    ∂ₜ(C) == D*Δ(C)
  end
  # Define Domain
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
  s,sd = circle(7, 500)
  # Create initial data.
  Csin = map(p -> sin(p[1]), point(s))
  Ccos = map(p -> cos(p[1]), point(s))
  C = stack([Csin, Ccos])
  u₀ = ComponentArray(C=Csin,)
  constants_and_parameters = (D = 0.001,)
  # Run
  function generate(sd, my_symbol; hodge=GeometricHodge()) end
  sim = eval(gensim(Heat,dimension=1))
  fₘ = sim(sd, nothing)
  tₑ = 1.15
  ode_prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
  ens_prob = EnsembleProblem(ode_prob,
    prob_func = (prob, i, repeat) ->
      remake(prob, u0=ComponentArray(C=C[:,i])))
  soln = solve(ens_prob, Tsit5(); trajectories=2)
  @test all(soln[1].u[1] .== Csin)
  @test all(soln[1].u[1] .!= Ccos)
  @test all(soln[2].u[1] .!= Csin)
  @test all(soln[2].u[1] .== Ccos)
end

@testset "Large Summations" begin
# Elementwise summations of more than 32 variables are not pre-compiled by our
# host language.

# Test that (.+)(...) is generated for small sums.
SmallSum = @decapode begin
  (A00, A01, A02, A03, A04, A05, A06, A07, A08, A09,
   A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,
   A20, A21, A22, A23, A24, A25, A26, A27, A28, A29,
   A30, A31, A32)::Form0

  ∂ₜ(A00) ==
          A01 + A02 + A03 + A04 + A05 + A06 + A07 + A08 + A09 +
    A10 + A11 + A12 + A13 + A14 + A15 + A16 + A17 + A18 + A19 +
    A20 + A21 + A22 + A23 + A24 + A25 + A26 + A27 + A28 + A29 +
    A30 + A31 + A32
end
needle = "A00̇ .= (.+)(A01, A02, A03, A04, A05, A06, A07, A08, A09, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19, A20, A21, A22, A23, A24, A25, A26, A27, A28, A29, A30, A31, A32)"
haystack = string(gensim(SmallSum))
@test occursin(needle, haystack)

# Test that sum([...]) is generated for large sums.
LargeSum = @decapode begin
  (A00, A01, A02, A03, A04, A05, A06, A07, A08, A09,
   A10, A11, A12, A13, A14, A15, A16, A17, A18, A19,
   A20, A21, A22, A23, A24, A25, A26, A27, A28, A29,
   A30, A31, A32, A33)::Form0

  ∂ₜ(A00) ==
          A01 + A02 + A03 + A04 + A05 + A06 + A07 + A08 + A09 +
    A10 + A11 + A12 + A13 + A14 + A15 + A16 + A17 + A18 + A19 +
    A20 + A21 + A22 + A23 + A24 + A25 + A26 + A27 + A28 + A29 +
    A30 + A31 + A32 + A33
end
needle = "A00̇ .= sum([A01, A02, A03, A04, A05, A06, A07, A08, A09, A10, A11, A12, A13, A14, A15, A16, A17, A18, A19, A20, A21, A22, A23, A24, A25, A26, A27, A28, A29, A30, A31, A32, A33])"
haystack = string(gensim(LargeSum))
@test occursin(needle, haystack)

end
