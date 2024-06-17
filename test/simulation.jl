using ACSets
using CombinatorialSpaces
using ComponentArrays
using Decapodes
using DiagrammaticEquations
using Distributions
using GeometryBasics: Point2, Point3
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using Test
Point3D = Point3{Float64}
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
@test Decapodes.get_vars_code(DiffusionWithConstant, [:k], Float64, CPUTarget()).args[2] == :(k = p.k)
@test infer_state_names(DiffusionWithConstant) == [:C, :k]

@test infer_state_names(DiffusionWithParameter) == [:C, :k]
@test Decapodes.get_vars_code(DiffusionWithParameter, [:k], Float64, CPUTarget()).args[2] == :(k = p.k(t))

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

f = evalsim(DiffusionWithLiteral, can_prealloc=false)
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
  infer_types!(budyko_sellers, op1_inf_rules_1D, op2_inf_rules_1D)
  resolve_overloads!(budyko_sellers, op1_res_rules_1D, op2_res_rules_1D)

  # This test ensures that the next one does not break, since it depends on
  # arbitrary internal variable naming.
  @test budyko_sellers[only(incident(budyko_sellers, Symbol("•1"), :name)), :type] == :DualForm0
  # A dual 0-form consists of ne(s) floats.
  @test occursin("var\"__•1\" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))",
    repr(gensim(budyko_sellers, dimension=1)))
end

@testset "Gensim Transformations" begin

  function checkForContractionInGensim(d::SummationDecapode)
    results = []
    block = gensim(d).args[2].args[2].args[5]
    for line in 2:length(block.args)
      push!(results, block.args[line].args[1])
    end

    return results
  end

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
    (A,C)::Form0
    (D)::Form2

    B == ∂ₜ(A)
    D == ∂ₜ(C)

    B == ⋆(⋆(A))
    D == d(d(C))
  end
  @test 4 == length(checkForContractionInGensim(single_contract))

  sim = eval(gensim(single_contract))
  f = sim(earth, default_dec_generate)
  A = 2 * ones(nv(earth))
  C = ones(nv(earth))
  u = ComponentArray(A=A, C=C)
  du = ComponentArray(A=zeros(nv(earth)), C=zeros(ntriangles(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A ≈ 2 * ones(nv(earth))
  @test du.C == zeros(ntriangles(earth))

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
  @test 4 == length(checkForContractionInGensim(single_contract))

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
  @test 4 == length(checkForContractionInGensim(single_contract))

  let sim = eval(gensim(contract_with_op2))
    f = sim(earth, default_dec_generate)
    A = 3 * ones(nv(earth))
    E_dec = ones(nv(earth))
    u = ComponentArray(A=A, E=E_dec)
    du = ComponentArray(A=zeros(ntriangles(earth)), E=zeros(nv(earth)))
    constants_and_parameters = ()
    f(du, u, constants_and_parameters, 0)

    @test du.A == zeros(ntriangles(earth))
    @test du.E ≈ 9 * ones(nv(earth))
  end

  let sim = eval(gensim(contract_with_op2, can_prealloc=false))
    f = sim(earth, default_dec_generate)
    A = 3 * ones(nv(earth))
    E_dec = ones(nv(earth))
    u = ComponentArray(A=A, E=E_dec)
    du = ComponentArray(A=zeros(ntriangles(earth)), E=zeros(nv(earth)))
    constants_and_parameters = ()
    f(du, u, constants_and_parameters, 0)

    @test du.A == zeros(ntriangles(earth))
    @test du.E ≈ 9 * ones(nv(earth))
  end


  # Testing contract lines beyond the initial value
  later_contraction = @decapode begin
    (A)::Form0

    D == ∂ₜ(A)

    B == A * A
    D == ⋆(⋆(B))
  end
  @test 4 == length(checkForContractionInGensim(single_contract))

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
  @test 0 == length(checkForContractionInGensim(no_contraction))

  sim = eval(gensim(no_contraction))
  f = sim(earth, default_dec_generate)
  A = [i for i in 1:nv(earth)]
  u = ComponentArray(A=A)
  du = ComponentArray(A=zeros(ne(earth)))
  constants_and_parameters = ()
  f(du, u, constants_and_parameters, 0)

  @test du.A == d(0, earth) * A

  # Testing no contraction of unallowed operators
  no_unallowed = @decapode begin
    (A)::Form0
    (D)::Form1

    D == ∂ₜ(A)
    D == d(k(A))
  end
  @test 0 == length(checkForContractionInGensim(no_unallowed))

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

  @test du.A == d(0, earth) * 20 * A

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

  let sim = eval(gensim(wedges01))
    f = sim(earth, default_dec_generate)
    A = ones(nv(earth))
    B = 2 * ones(nv(earth))
    C = 3 * ones(ne(earth))
    u = ComponentArray(A=A, B=B, C=C)
    du = ComponentArray(A=zeros(ne(earth)), B=zeros(ne(earth)), C=zeros(ne(earth)))

    constants_and_parameters = ()
    f(du, u, constants_and_parameters, 0)

    @test du.A == du.B == du.C
  end

  let sim = eval(gensim(wedges01, can_prealloc=false))
    f = sim(earth, default_dec_generate)
    A = ones(nv(earth))
    B = 2 * ones(nv(earth))
    C = 3 * ones(ne(earth))
    u = ComponentArray(A=A, B=B, C=C)
    du = ComponentArray(A=zeros(ne(earth)), B=zeros(ne(earth)), C=zeros(ne(earth)))

    constants_and_parameters = ()
    f(du, u, constants_and_parameters, 0)

    @test du.A == du.B == du.C
  end


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
      B == ⋆(A)
    end
    g = gensim(DiagonalInvHodge1)
    @test gensim(DiagonalInvHodge1).args[2].args[2].args[3].args[2].args[2].args[3].value == :⋆₁⁻¹
    sim = eval(g)

    # Test that no error is thrown here
    f = sim(line, default_dec_generate, DiagonalHodge())
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

end

@testset "Allocations" begin
# Test the heat equation Decapode has expected memory allocation.
Heat = @decapode begin
  C::Form0
  D::Constant
  ∂ₜ(C) == D*Δ(C)
end
sim = eval(gensim(Heat))
s = loadmesh(Icosphere(1))
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())
f = sim(sd,nothing)
u₀ = ComponentArray(
  C = map(x -> x[3], point(sd)))
p = (D=1e-1,)
du = copy(u₀)
# The first call to the function makes many allocations.
_ = @allocations f(du, u₀, p, (0,1.0)) # 55259
_ = @allocated f(du, u₀, p, (0,1.0)) # 3962696
# Test that subsequent calls make a reasonable amount.
bytes = @allocated f(du, u₀, p, (0,1.0))
nallocs = @allocations f(du, u₀, p, (0,1.0))
@test nallocs == 3
@test bytes == 80
end
