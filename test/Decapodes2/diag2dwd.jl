using Test
using Catlab
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Present
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using LinearAlgebra
using MLStyle
using Base.Iterators

using Decapodes
import Decapodes: DecaExpr

# @present DiffusionSpace2D(FreeExtCalc2D) begin
#   X::Space
#   k::Hom(Form1(X), Form1(X)) # diffusivity of space, usually constant (scalar multiplication)
#   proj₁_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
#   proj₂_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
#   sum₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
#   prod₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
# end


# Diffusion = @decapode DiffusionSpace2D begin
#     (C, Ċ₁, Ċ₂)::Form0{X}
#     Ċ₁ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(k(d₀{X}(C)))))
#     Ċ₂ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(d₀{X}(C))))
#     ∂ₜ{Form0{X}}(C) == Ċ₁ + Ċ₂
# end

# Tests
#######

# Construct roughly what the @decapode macro should return for Diffusion
js = [Judge(Var(:C), :Form0, :X), 
      Judge(Var(:Ċ₁), :Form0, :X),
      Judge(Var(:Ċ₂), :Form0, :X)
]
# TODO: Do we need to handle the fact that all the functions are parameterized by a space?
eqs = [Eq(Var(:Ċ₁), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Var(:Ċ₂), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus([Var(:Ċ₁), Var(:Ċ₂)]))
]
diffusion_d = DecaExpr(js, eqs)
# diffusion_cset = Decapode(diffusion_d)
diffusion_cset_named = SummationDecapode(diffusion_d)
# A test with expressions on LHS (i.e. temporary variables must be made)
# TODO: we should be intelligent and realize that the LHS of the first two
# equations are the same and so can share a new variable
eqs = [Eq(Plus([Var(:Ċ₁), Var(:Ċ₂)]), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Plus([Var(:Ċ₁), Var(:Ċ₂)]), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus([Var(:Ċ₁), Var(:Ċ₂)]))    
]
test_d = DecaExpr(js, eqs)
# test_cset = Decapode(test_d)
test_cset_named = SummationDecapode(test_d)

# TODO: Write tests for recursive expressions

all(isassigned(test_cset_named[:name], i) for i in parts(test_cset_named,:Var))

sup_js = js = [Judge(Var(:C), :Form0, :X), 
Judge(Var(:ϕ₁), :Form0, :X),
Judge(Var(:ϕ₂), :Form0, :X)
]
sup_eqs = [Eq(Var(:ϕ₁), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Var(:ϕ₂), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus([Var(:ϕ₁), Var(:ϕ₂)]))
]
sup_d = DecaExpr(sup_js, sup_eqs)
# sup_cset = Decapode(sup_d)
sup_cset_named = SummationDecapode(sup_d)


Decapodes.compile(diffusion_cset_named, [:C,])
Decapodes.compile(test_cset_named, [:C,])
Decapodes.compile(sup_cset_named, [:C,])

term(:(∧₀₁(C,V)))

@testset "Term Construction" begin
    @test term(:(Ċ)) == Var(:Ċ)
    @test_throws ErrorException term(:(∂ₜ{Form0}))
    @test term(Expr(:ϕ)) == Var(:ϕ)
    @test typeof(term(:(d₀(C)))) == App1
    @test typeof(term(:(∘(k, d₀)(C)))) == AppCirc1
    # @test term(:(∘(k, d₀)(C))) == AppCirc1([:k, :d₀], Var(:C)) #(:App1, ((:Circ, :k, :d₀), Var(:C)))
    # @test term(:(∘(k, d₀{X})(C))) == (:App1, ((:Circ, :k, :(d₀{X})), Var(:C)))
    @test_throws MethodError term(:(Ċ == ∘(⋆₀⁻¹{X}, dual_d₁{X}, ⋆₁{X})(ϕ)))
    @test term(:(∂ₜ(C))) == Tan(Var(:C))
    # @test term(:(∂ₜ{Form0}(C))) == App1(:Tan, Var(:C))
end

@testset "Recursive Expr" begin
  Recursion = quote
    x::Form0{X}
    y::Form0{X}
    z::Form0{X}

    ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
    y == F(f2(x), ρ(x,z))
  end

  recExpr = parse_decapode(Recursion)
  rdp = SummationDecapode(recExpr)
  show(rdp)

  @test nparts(rdp, :Var) == 9
  @test nparts(rdp, :TVar) == 1
  @test nparts(rdp, :Op1) == 5
  @test nparts(rdp, :Op2) == 2
  @test nparts(rdp, :Σ) == 1
end
Recursion = quote
  x::Form0{X}
  y::Form0{X}
  z::Form0{X}

  ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
  y == F(f2(x), ρ(x,z))
end

recExpr = parse_decapode(Recursion)
rdp = SummationDecapode(recExpr)

@testset "Diffusion Diagram" begin
    DiffusionExprBody =  quote
        (C, Ċ)::Form0{X}
        ϕ::Form1{X}
    
        # Fick's first law
        ϕ ==  ∘(k, d₀)(C)
        # Diffusion equation
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    diffExpr = parse_decapode(DiffusionExprBody)
    ddp = SummationDecapode(diffExpr)
    to_graphviz(ddp)

    @test nparts(ddp, :Var) == 3
    @test nparts(ddp, :TVar) == 1
    @test nparts(ddp, :Op1) == 3
    @test nparts(ddp, :Op2) == 0
end


@testset "Advection Diagram" begin
    Advection = quote
        C::Form0{X}
        (V, ϕ)::Form1{X}

        ϕ == ∧₀₁(C,V)
    end

    advdecexpr = parse_decapode(Advection)
    advdp = SummationDecapode(advdecexpr)
    @test nparts(advdp, :Var) == 3
    @test nparts(advdp, :TVar) == 0
    @test nparts(advdp, :Op1) == 0
    @test nparts(advdp, :Op2) == 1
end

@testset "Superposition Diagram" begin
    Superposition = quote
        (C, Ċ)::Form0{X}
        (ϕ, ϕ₁, ϕ₂)::Form1{X}

        ϕ == ϕ₁ + ϕ₂
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    superexp = parse_decapode(Superposition)
    supdp = SummationDecapode(superexp)
    @test nparts(supdp, :Var) == 5
    @test nparts(supdp, :TVar) == 1
    @test nparts(supdp, :Op1) == 2
    @test nparts(supdp, :Op2) == 0
    @test nparts(supdp, :Σ) == 1
    @test nparts(supdp, :Summand) == 2
end

@testset "AdvectionDiffusion Diagram" begin
    AdvDiff = quote
        (C, Ċ)::Form0{X}
        (V, ϕ, ϕ₁, ϕ₂)::Form1{X}

        # Fick's first law
        ϕ₁ ==  (k ∘ d₀)(C)
        # Diffusion equation
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ϕ₂ == ∧₀₁(C,V)

        ϕ == ϕ₁ + ϕ₂
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    advdiff = parse_decapode(AdvDiff)
    advdiffdp = SummationDecapode(advdiff)
    @test nparts(advdiffdp, :Var) == 6
    @test nparts(advdiffdp, :TVar) == 1
    @test nparts(advdiffdp, :Op1) == 4
    @test nparts(advdiffdp, :Op2) == 1
    @test nparts(advdiffdp, :Σ) == 1
    @test nparts(advdiffdp, :Summand) == 2
end

@testset "Type Inference" begin
  # The type of the tgt of ∂ₜ is inferred.
  Test1 = quote
    C::Form0{X}
    ∂ₜ(C) == C
  end
  t1 = SummationDecapode(parse_decapode(Test1))
  infer_types!(t1)

  # We use set equality because we do not care about the order of the Var table.
  names_types_1 = Set(zip(t1[:name], t1[:type]))
  names_types_expected_1 = Set([(:Ċ, :Form0), (:C, :Form0)])
  @test issetequal(names_types_1, names_types_expected_1)

  # The type of the src of ∂ₜ is inferred.
  Test2 = quote
    C::infer{X}
    ∂ₜ(C) == C
  end
  t2 = SummationDecapode(parse_decapode(Test2))
  t2[:type][only(incident(t2, :Ċ, :name))] = :Form0
  infer_types!(t2)

  names_types_2 = Set(zip(t2[:name], t2[:type]))
  names_types_expected_2 = Set([(:Ċ, :Form0), (:C, :Form0)])
  @test issetequal(names_types_2, names_types_expected_2)

  # The type of the tgt of d is inferred.
  Test3 = quote
    C::Form0{X}
    D::infer{X}
    E::infer{X}
    #C::infer{X}
    #D::Form1{X}
    #E::infer{X}
    D == d(C)
    E == d(D)
  end
  t3 = SummationDecapode(parse_decapode(Test3))
  #t3_inferred = infer_types!(t3)
  infer_types!(t3)

  names_types_3 = Set(zip(t3[:name], t3[:type]))
  names_types_expected_3 = Set([(:C, :Form0), (:D, :Form1), (:E, :Form2)])
  @test issetequal(names_types_3, names_types_expected_3)

  # The type of the src and tgt of d is inferred.
  Test4 = quote
    C::infer{X}
    D::Form1{X}
    E::infer{X}
    D == d(C)
    E == d(D)
  end
  t4 = SummationDecapode(parse_decapode(Test4))
  #t4_inferred = infer_types!(t4)
  infer_types!(t4)

  names_types_4 = Set(zip(t4[:name], t4[:type]))
  names_types_expected_4 = Set([(:C, :Form0), (:D, :Form1), (:E, :Form2)])
  @test issetequal(names_types_4, names_types_expected_4)

  # The type of the src of d is inferred.
  Test5 = quote
    C::infer{X}
    D::Form1{X}
    D == d(C)
  end
  t5 = SummationDecapode(parse_decapode(Test5))
  infer_types!(t5)

  names_types_5 = Set(zip(t5[:name], t5[:type]))
  names_types_expected_5 = Set([(:C, :Form0), (:D, :Form1)])
  @test issetequal(names_types_5, names_types_expected_5)

  # The type of the src of d is inferred.
  Test6 = quote
    A::Form0{X}
    (B,C,D,E,F)::infer{X}
    B == d(A)
    C == d(B)
    D == ⋆(C)
    E == d(D)
    F == d(E)
  end
  t6 = SummationDecapode(parse_decapode(Test6))
  infer_types!(t6)

  names_types_6 = Set(zip(t6[:name], t6[:type]))
  names_types_expected_6 = Set([
    (:A, :Form0),     (:B, :Form1),     (:C, :Form2),
    (:F, :DualForm2), (:E, :DualForm1), (:D, :DualForm0)])
  @test issetequal(names_types_6, names_types_expected_6)

  # The type of a summand is inferred.
  Test7 = quote
    A::Form0{X}
    (B,C,D,E,F)::infer{X}
    A == B+C+D+E+F
  end
  t7 = SummationDecapode(parse_decapode(Test7))
  infer_types!(t7)

  types_7 = Set(t7[:type])
  types_expected_7 = Set([:Form0])
  @test issetequal(types_7, types_expected_7)

  # The type of a sum is inferred.
  Test8 = quote
    B::Form0{X}
    (A,C,D,E,F)::infer{X}
    A == B+C+D+E+F
  end
  t8 = SummationDecapode(parse_decapode(Test8))
  infer_types!(t8)

  types_8 = Set(t8[:type])
  types_expected_8 = Set([:Form0])
  @test issetequal(types_8, types_expected_8)

  # Type inference of op1 propagates through sums.
  Test8 = quote
    Γ::Form0{X}
    (A,B,C,D,E,F,Θ)::infer{X}
    B == d(Γ)
    A == B+C+D+E+F
    Θ == d(A)
  end
  t8 = SummationDecapode(parse_decapode(Test8))
  infer_types!(t8)

  names_types_8 = Set(zip(t8[:name], t8[:type]))
  names_types_expected_8 = Set([
    (:Γ, :Form0),
    (:A, :Form1), (:B, :Form1), (:C, :Form1), (:D, :Form1), (:E, :Form1), (:F, :Form1),
    (:Θ, :Form2)])
  @test issetequal(names_types_8, names_types_expected_8)

  function makeInferPathDeca(log_cycles; infer_path = false)
    cycles = 2 ^ log_cycles
    num_nodes = 6 * cycles
    num_ops = num_nodes - 1

    if(infer_path)
      type_arr = vcat(:Form0, fill(:infer, num_nodes - 1))
    else
      type_arr = repeat([:Form0, :Form1, :Form2, :DualForm0, :DualForm1, :DualForm2], cycles)
    end

    InferPathTest = @acset SummationDecapode{Any, Any, Symbol} begin
      Var = num_nodes
      type = type_arr
      name = map(x -> Symbol("•$x"), 1:num_nodes)

      Op1 = num_ops
      src = 1:num_ops
      tgt = 2:num_ops + 1
      op1 = repeat([:d, :d, :⋆, :d, :d, :⋆], cycles)[1:num_ops]
    end
  end
  # log_cycles > 6 starts to be slow
  log_cycles = 3

  decapode_to_infer = makeInferPathDeca(log_cycles, infer_path = true)
  decapode_expected = makeInferPathDeca(log_cycles)
  infer_types!(decapode_to_infer)
  @test decapode_to_infer == decapode_expected

  # Inference can happen "through" composed ops.
  Test9 = quote
    A::Form0{X}
    B::infer{X}
    B == ∘(d,d,⋆,d,d)(A)
  end
  t9 = SummationDecapode(parse_decapode(Test9))
  expand_compositions!(t9)
  infer_types!(t9)

  names_types_9 = Set(zip(t9[:name], t9[:type]))
  names_types_expected_9 = Set([
    (:A, :Form0),     (Symbol('•', 1), :Form1),     (Symbol('•', 2), :Form2),
    (:B, :DualForm2), (Symbol('•', 4), :DualForm1), (Symbol('•', 3), :DualForm0)])
  @test issetequal(names_types_9, names_types_expected_9)
end

@testset "Overloading Resolution" begin
  # d overloading is resolved.
  Test1 = quote
    A::Form0{X}
    B::Form1{X}
    C::Form2{X}
    D::DualForm0{X}
    E::DualForm1{X}
    F::DualForm2{X}
    B == d(A)
    C == d(B)
    E == d(D)
    F == d(E)
  end
  t1 = SummationDecapode(parse_decapode(Test1))
  resolve_overloads!(t1)

  op1s_1 = t1[:op1]
  op1s_expected_1 = [:d₀ , :d₁, :dual_d₀, :dual_d₁]
  @test op1s_1 == op1s_expected_1
  
  # ⋆ overloading is resolved.
  Test2 = quote
    C::Form0{X}
    D::DualForm2{X}
    E::Form1{X}
    F::DualForm1{X}
    G::Form2{X}
    H::DualForm0{X}
    D == ⋆(C)
    C == ⋆(D)
    E == ⋆(F)
    F == ⋆(E)
    G == ⋆(H)
    H == ⋆(G)
  end
  t2 = SummationDecapode(parse_decapode(Test2))
  resolve_overloads!(t2)

  op1s_2 = t2[:op1]
  # Note: The Op1 table of the decapode created does not have the functions
  # listed in the order in which they appear in Test2.
  op1s_expected_2 = [:⋆₀ , :⋆₀⁻¹ , :⋆₁⁻¹ , :⋆₁ , :⋆₂⁻¹ , :⋆₂]
  @test op1s_2 == op1s_expected_2
  
  # All overloading on the de Rahm complex is resolved.
  Test3 = quote
    A::Form0{X}
    B::Form1{X}
    C::Form2{X}
    D::DualForm0{X}
    E::DualForm1{X}
    F::DualForm2{X}
    B == d(A)
    C == d(B)
    E == d(D)
    F == d(E)
    F == ⋆(A)
    A == ⋆(F)
    E == ⋆(B)
    B == ⋆(E)
    D == ⋆(C)
    C == ⋆(D)
  end
  t3 = SummationDecapode(parse_decapode(Test3))
  resolve_overloads!(t3)

  op1s_3 = t3[:op1]
  # Note: The Op1 table of the decapode created does not have the functions
  # listed in the order in which they appear in Test2.
  op1s_expected_3 = [:d₀ , :d₁, :dual_d₀, :dual_d₁, :⋆₀ , :⋆₀⁻¹ , :⋆₁ , :⋆₁⁻¹ , :⋆₂ , :⋆₂⁻¹]
  @test op1s_3 == op1s_expected_3
end

