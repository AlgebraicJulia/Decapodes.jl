using Decapodes
import Decapodes: compile, gensim, infer_states, infer_state_names

using Catlab
using Catlab.CategoricalAlgebra

using Test

using MLStyle
using CombinatorialSpaces
using GeometryBasics: Point3
using LinearAlgebra
using Distributions
using ComponentArrays
using OrdinaryDiffEq

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

DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ == ∘(k, d₀)(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
add_constant!(ddp, :k)
@test nparts(ddp, :Var) == 4

DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}
    k::Constant{ℝ}


    # Fick's first law
    ϕ == k * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)

dec_matrices = Vector{Symbol}()
alloc_vectors = Vector{Decapodes.AllocVecCall}()

@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k)
@test infer_state_names(ddp) == [:C, :k]

torus = loadmesh(Torus_30x10())
c_dist = MvNormal([5, 5], LinearAlgebra.Diagonal(map(abs2, [1.5, 1.5])))
c = [pdf(c_dist, [p[1], p[2]]) for p in torus[:point]]

u₀ = ComponentArray(C=c)
du = ComponentArray(C=zero(c))

f = eval(gensim(expand_operators(ddp)))
fₘₛ = f(torus, generate)

DiffusionExprBody =  quote
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}
  k::Parameter{ℝ}

  # Fick's first law
  ϕ == k * d₀(C)
  # Diffusion equation
  Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
  ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
dec_matrices2 = Vector{Symbol}()
alloc_vectors2 = Vector{Decapodes.AllocVecCall}()

@test infer_state_names(ddp) == [:C, :k]
@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k(t))

f = eval(gensim(expand_operators(ddp)))
fₘₚ = f(torus, generate)

@test norm(fₘₛ(du, u₀, (k=2.0,), 0)  - fₘₚ(du, u₀, (k=t->2.0,), 0)) < 1e-4

DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ == 3 * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)

@test infer_state_names(ddp) == [:C]
# TODO: Fix proper Expr equality, the Float64 does not equate here
# @test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2] == :(var"3"::Float64 = 3.0)
@test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2].args[1] == :(var"3")
@test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2].args[2].args[3] == :(3.0)


f = eval(gensim(expand_operators(ddp)))
fₘₚ = f(torus, generate)

@test norm(fₘₛ(du, u₀, (k=2.0,), 0)  - fₘₚ(du, u₀, (k=t->2.0,), 0)) < 1e-4

DiffusionExprBody =  quote
    (C, Ċ)::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ == 3 * d₀(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)

@test infer_state_names(ddp) == [:C]
# TODO: Fix proper Expr equality, the Float64 does not equate here
# @test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2] == :(var"3"::Float64 = 3.0)
@test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2].args[1] == :(var"3")
@test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2].args[2].args[3] == :(3.0)

f = eval(gensim(expand_operators(ddp)))
fₘₚ = f(torus, generate)

@test norm(fₘₛ(du, u₀, (k=2.0,), 0)  - fₘₚ(du, u₀, (k=t->2.0,), 0)) < 1e-4
 
# to solve the ODE over a duration, use the ODEProblem from OrdinaryDiffEq
# tₑ = 10
# using OrdinaryDiffEq
# prob = ODEProblem(fₘₛ,u₀, (0,tₑ), (k=2.0,))
# soln = solve(prob, Tsit5())
end

# Testing done based on the original gensim
Point3D = Point3{Float64}
flatten(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

# Testing Advection-Diffusion
#= @testset "Advection-Diffusion Simulation" begin

  function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :k => x->2000x
      :d₀ => test_differential(0, sd)
      :⋆₁ => test_hodge(1, sd, hodge)
      :⋆₀⁻¹ => test_inverse_hodge(0, sd, hodge)
      :dual_d₁ => test_dual_differential(1, sd)  
      :(-) => x-> -x
      :plus => (+)
      _ => default_dec_generate_2D(sd, my_symbol, hodge)
    end

    return (args...) ->  op(args...)
  end

  
  RADIUS = 6371+90
  primal_earth = loadmesh(Icosphere(1, RADIUS))
  nploc = argmax(x -> x[3], primal_earth[:point])
  primal_earth[:edge_orientation] = false
  orient!(primal_earth)
  earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
  subdivide_duals!(earth, Circumcenter())
  
  AdvDiff = quote
      C::Form0{X}
      Ċ::Form0{X}
      V::Form1{X}
      ϕ::Form1{X}
      ϕ₁::Form1{X}
      ϕ₂::Form1{X}
      # starC::DualForm2{X}
      # lvc::Form1{X}
      # Fick's first law
      ϕ₁ ==  ∘(d₀,k,⋆₁)(C)
      # Advective Flux
      ϕ₂ == -(L₀(V, C))
      # Superposition Principle
      ϕ == ϕ₁ + ϕ₂
      # Conservation of Mass
      Ċ == ∘(dual_d₁,⋆₀⁻¹)(ϕ)
      ∂ₜ(C) == Ċ
  end

  advdiff = parse_decapode(AdvDiff)
  advdiffdp = SummationDecapode(advdiff)

  function old_simulate(mesh, operators)
    begin
        d₀ = generate(mesh, :d₀)
        k = generate(mesh, :k)
        (⋆₁) = generate(mesh, :⋆₁)
        (⋆₀) = generate(mesh, :⋆₀)
        (-) = generate(mesh, :-)
        dual_d₁ = generate(mesh, :dual_d₁)
        (⋆₀⁻¹) = generate(mesh, :⋆₀⁻¹)
        L₀ = operators(mesh, :L₀)
        plus = operators(mesh, :plus)
    end
    return begin
            f(du, u, p, t) = begin
                    begin
                        C = (u.C)
                        V = (u.V)
                    end
                    ϕ₁ = (∘(⋆₁, k, d₀))(C)
                    var"•1" = L₀(V, C)
                    ϕ₂ = -var"•1"
                    ϕ = plus(ϕ₁, ϕ₂)
                    Ċ = ((⋆₀⁻¹) ∘ dual_d₁)(ϕ)
                    du .= 0.0
                    begin
                        (du.C) .= Ċ
                    end
                end
        end
  end

  fₙ = old_simulate(earth, generate)

  new_sim = evalsim(advdiffdp, [:C, :V])
  fₘ = new_sim(earth, generate)

  # Running with the old gensim simulation
  velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4))
  begin
    c_dist = MvNormal(nploc[[1,2]], 100[1, 1])
    c = [pdf(c_dist, [p[1], p[2]]./√RADIUS) for p in earth[:point]]

    vmag = 500
    v = flatten(velocity, earth)

    u₀ = ComponentArray(C=c,V=collect(v))
    tₑ = 6
    prob = ODEProblem(fₙ,u₀,(0, tₑ))
    old_soln = solve(prob, Tsit5())
  end

  # Running with the new gensim simulation
  begin
    c_dist = MvNormal(nploc[[1,2]], 100[1, 1])
    c = [pdf(c_dist, [p[1], p[2]]./√RADIUS) for p in earth[:point]]

    vmag = 500
    v = flatten(velocity, earth)

    u₀ = ComponentArray(C=c,V=collect(v))
    tₑ = 6
    prob = ODEProblem(fₘ,u₀,(0, tₑ))
    new_soln = solve(prob, Tsit5())
  end

  @test old_soln.u ≈ new_soln.u
end =#

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
  result = f(u, u, constants_and_parameters, 0)

  @test result == -1 * ones(nv(earth))


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
  
  sim = eval(gensim(contract_with_op2))
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
  
  sim = eval(gensim(wedges01))

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

  @test all(isapprox.(du.A, zeros(ntriangles(earth)); atol = 1e-300))
  @test all(isapprox.(du.B, zeros(ntriangles(earth)); atol = 1e-300))
  @test all(isapprox.(du.A, du.B; atol = 1e-300))

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