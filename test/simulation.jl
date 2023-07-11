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
using MultiScaleArrays
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
    _ => default_dec_generate(sd, my_symbol)
    # :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length(args[1]))"); op(args...);end
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
compile(expand_operators(ddp), [:C, :k], dec_matrices, alloc_vectors)

@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k)
@test infer_state_names(ddp) == [:C, :k]

gensim(ddp)

torus = loadmesh(Torus_30x10())
c_dist = MvNormal([5, 5], LinearAlgebra.Diagonal(map(abs2, [1.5, 1.5])))
c = [pdf(c_dist, [p[1], p[2]]) for p in torus[:point]]

u₀ = construct(PhysicsState, [VectorForm(c)],Float64[], [:C])
du = construct(PhysicsState, [VectorForm(zero(c))],Float64[], [:C])

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

compile(expand_operators(ddp), [:C, :k], dec_matrices2, alloc_vectors2)


@test infer_state_names(ddp) == [:C, :k]
@test Decapodes.get_vars_code(ddp, [:k]).args[2] == :(k = p.k(t))
gensim(ddp)

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
gensim(ddp)

@test infer_state_names(ddp) == [:C]
@test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2] == :(var"3" = 3.0)

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
gensim(ddp)

@test infer_state_names(ddp) == [:C]
@test Decapodes.get_vars_code(ddp, [Symbol("3")]).args[2] == :(var"3" = 3.0)

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
@testset "Advection-Diffusion Simulation" begin

  function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :k => x->2000x
      :d₀ => test_differential(0, sd)
      :⋆₁ => test_hodge(1, sd, hodge)
      :⋆₀⁻¹ => test_inverse_hodge(0, sd, hodge)
      :dual_d₁ => test_dual_differential(1, sd)  
      :(-) => x-> -x
      :plus => (+)
      _ => default_dec_generate(sd, my_symbol, hodge)
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
      ϕ == plus(ϕ₁ , ϕ₂)
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
        (-) = generate(mesh, :-)
        dual_d₁ = generate(mesh, :dual_d₁)
        (⋆₀⁻¹) = generate(mesh, :⋆₀⁻¹)
        L₀ = operators(mesh, :L₀)
        plus = operators(mesh, :plus)
    end
    return begin
            f(du, u, p, t) = begin
                    begin
                        C = (findnode(u, :C)).values
                        V = (findnode(u, :V)).values
                    end
                    ϕ₁ = (∘(⋆₁, k, d₀))(C)
                    var"•1" = L₀(V, C)
                    ϕ₂ = -var"•1"
                    ϕ = plus(ϕ₁, ϕ₂)
                    Ċ = ((⋆₀⁻¹) ∘ dual_d₁)(ϕ)
                    du .= 0.0
                    begin
                        (findnode(du, :C)).values .= Ċ
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

    u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(collect(v))],Float64[], [:C, :V])
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

    u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(collect(v))],Float64[], [:C, :V])
    tₑ = 6
    prob = ODEProblem(fₘ,u₀,(0, tₑ))
    new_soln = solve(prob, Tsit5())
  end

  @test old_soln.u ≈ new_soln.u
end

# Testing Brusselator
@testset "Brusselator Simulation" begin

  function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :Δ₀ => test_laplace_de_rham(0, sd)
        _ => default_dec_generate(sd, my_symbol, hodge)
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
                        U = (findnode(u, :U)).values
                        V = (findnode(u, :V)).values
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
                        (findnode(du, :U)).values .= U̇
                        (findnode(du, :V)).values .= V̇
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

    u₀ = construct(PhysicsState, [VectorForm(U), VectorForm(V)],Float64[], [:U, :V])
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

      u₀ = construct(PhysicsState, [VectorForm(U), VectorForm(V)],Float64[], [:U, :V])
      tₑ = 11.5
      prob = ODEProblem(fₘ,u₀,(0, tₑ), constants_and_parameters)
      new_soln = solve(prob, Tsit5())
  end

  @test old_soln.u ≈ new_soln.u
end
