using Test
using Catlab
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Present
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.Graphics
using Catlab.Programs
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
  # Warning, this testing depends on the fact that varname, form information is 
  # unique within a decapode even though this is not enforced
  function get_name_type_pair(d::SummationDecapode)
    Set(zip(d[:name], d[:type]))
  end

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
  t9 = expand_operators(t9)
  infer_types!(t9)

  names_types_9 = Set(zip(t9[:name], t9[:type]))
  names_types_expected_9 = Set([
    (:A, :Form0),     (Symbol("•_1_", 1), :Form1),     (Symbol("•_1_", 2), :Form2),
    (:B, :DualForm2), (Symbol("•_1_", 4), :DualForm1), (Symbol("•_1_", 3), :DualForm0)])
  @test issetequal(names_types_9, names_types_expected_9)

  # Basic op2 inference using ∧
  Test10 = quote
    (A, B)::Form0{X}
    C::infer{X}

    D::Form0{X}
    E::Form1{X}
    (F, H)::infer{X}

    C == ∧(A, B)

    F == ∧(D, E)
    H == ∧(E, D)
  end
  t10 = SummationDecapode(parse_decapode(Test10))
  infer_types!(t10)

  names_types_10 = get_name_type_pair(t10)
  names_types_expected_10 = Set([(:F, :Form1), (:B, :Form0), (:C, :Form0), (:H, :Form1), (:A, :Form0), (:E, :Form1), (:D, :Form0)])
  @test issetequal(names_types_10, names_types_expected_10)

  # Basic op2 inference using L
  Test11 = quote
    (A, E)::Form1{X}
    B::Form0{X}
    (C, D)::infer{X}

    C == L(A, B)
    D == L(A, E)
  end
  t11 = SummationDecapode(parse_decapode(Test11))
  infer_types!(t11)

  names_types_11 = get_name_type_pair(t11)
  names_types_expected_11 = Set([(:A, :Form1), (:B, :Form0), (:C, :Form0), (:E, :Form1), (:D, :Form1)])
  @test issetequal(names_types_11, names_types_expected_11)

  # Basic op2 inference using i
  Test12 = quote
    (A, B)::Form1{X}
    C::infer{X}

    C == i(A, B)
  end
  t12 = SummationDecapode(parse_decapode(Test12))
  infer_types!(t12)

  names_types_12 = get_name_type_pair(t12)
  names_types_expected_12 = Set([(:A, :Form1), (:B, :Form1), (:C, :Form0)])
  @test issetequal(names_types_12, names_types_expected_12)

  #2D op2 inference using ∧
  Test13 = quote
    (A, B)::Form1{X}
    C::infer{X}

    D::Form0{X}
    E::Form2{X}
    (F, H)::infer{X}

    C == ∧(A, B)

    F == ∧(D, E)
    H == ∧(E, D)
  end
  t13 = SummationDecapode(parse_decapode(Test13))
  infer_types!(t13)

  names_types_13 = get_name_type_pair(t13)
  names_types_expected_13 = Set([(:E, :Form2), (:B, :Form1), (:C, :Form2), (:A, :Form1), (:F, :Form2), (:H, :Form2), (:D, :Form0)])
  @test issetequal(names_types_13, names_types_expected_13)

  # 2D op2 inference using L
  Test14 = quote
    A::Form1{X}
    B::Form2{X}
    C::infer{X}

    C == L(A, B)
  end
  t14 = SummationDecapode(parse_decapode(Test14))
  infer_types!(t14)

  names_types_14 = get_name_type_pair(t14)
  names_types_expected_14 = Set([(:C, :Form2), (:A, :Form1), (:B, :Form2)])
  @test issetequal(names_types_14, names_types_expected_14)

  # 2D op2 inference using i
  Test15 = quote
    A::Form1{X}
    B::Form2{X}
    C::infer{X}

    C == i(A, B)
  end
  t15 = SummationDecapode(parse_decapode(Test15))
  infer_types!(t15)

  names_types_15 = get_name_type_pair(t15)
  names_types_expected_15 = Set([(:A, :Form1), (:B, :Form2), (:C, :Form1)])
  @test issetequal(names_types_15, names_types_expected_15)

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

  Test4 = quote
    (A, C, F, H)::Form0{X}
    (B, D, E, G)::Form1{X}


    C == ∧(A, A)
    D == ∧(A, B)
    E == ∧(B, A)
    F == L(B, A)
    G == L(B, B)
    H == i(B, B)
  end
  t4 = SummationDecapode(parse_decapode(Test4))
  resolve_overloads!(t4)
  op2s_4 = t4[:op2]
  op2s_expected_4 = [:∧₀₀ , :∧₀₁, :∧₁₀, :L₀, :L₁, :i₁]
  @test op2s_4 == op2s_expected_4

  Test5 = quote
    A::Form0{X}
    (B, H)::Form1{X}
    (C, D, E, F, G)::Form2{X}

    D == ∧(B, B)
    E == ∧(A, C)
    F == ∧(C, A)
    G == L(B, C)
    H == i(B, C)
  end
  t5 = SummationDecapode(parse_decapode(Test5))
  resolve_overloads!(t5)
  op2s_5 = t5[:op2]
  op2s_expected_5 = [:∧₁₁, :∧₀₂, :∧₂₀, :L₂, :i₂]
  @test op2s_5 == op2s_expected_5
end

@testset "Type Inference and Overloading Resolution Integration" begin
  # Momentum-formulation of Navier Stokes on sphere
  DiffusionExprBody = quote
    (T, Ṫ)::Form0{X}
    ϕ::DualForm1{X}
    k::Constant{X}
    # Fick's first law
    ϕ ==  ⋆(k*d(T))
    # Diffusion equation
    Ṫ ==  ⋆(d(ϕ))
  end
  Diffusion = SummationDecapode(parse_decapode(DiffusionExprBody))
  AdvectionExprBody = quote
    (M,V)::Form1{X}  #  M = ρV
    (ρ, p, T, Ṫ)::Form0{X}
    V == M/avg(ρ)
    ρ == p / R₀(T)
    Ṫ == neg(⋆(L(V, ⋆(T))))
  end
  Advection = SummationDecapode(parse_decapode(AdvectionExprBody))
  SuperpositionExprBody = quote
    (T, Ṫ, Ṫ₁, Ṫₐ)::Form0{X}
    Ṫ == Ṫ₁ + Ṫₐ
    ∂ₜ(T) == Ṫ 
  end
  Superposition = SummationDecapode(parse_decapode(SuperpositionExprBody))
  compose_continuity = @relation () begin
    diffusion(T, Ṫ₁)
    advection(M, ρ, P, T, Ṫₐ)
    superposition(T, Ṫ, Ṫ₁, Ṫₐ)
  end
  continuity_cospan = oapply(compose_continuity,
                  [Open(Diffusion, [:T, :Ṫ]),
                   Open(Advection, [:M, :ρ, :p, :T, :Ṫ]),
                   Open(Superposition, [:T, :Ṫ, :Ṫ₁, :Ṫₐ])])

  continuity = apex(continuity_cospan)
  NavierStokesExprBody = quote
    (M, Ṁ, G, V)::Form1{X}
    (T, ρ, p, ṗ)::Form0{X}
    (two,three,kᵥ)::Constant{X}
    V == M/avg(ρ)
    Ṁ == neg(L(V, V))*avg(ρ) + 
          kᵥ*(Δ(V) + d(δ(V))/three) +
          d(i(V, V)/two)*avg(ρ) +
          neg(d(p)) +
          G*avg(ρ)
    ∂ₜ(M) == Ṁ
    ṗ == neg(⋆(L(V, ⋆(p)))) # *Lie(3Form) = Div(*3Form x v) --> conservation of pressure
    ∂ₜ(p) == ṗ
  end
  NavierStokes = SummationDecapode(parse_decapode(NavierStokesExprBody))
  compose_heatXfer = @relation () begin
    continuity(M, ρ, P, T)
    navierstokes(M, ρ, P, T)
  end
  heatXfer_cospan = oapply(compose_heatXfer,
                  [Open(continuity, [:M, :ρ, :P, :T]),
                   Open(NavierStokes, [:M, :ρ, :p, :T])])
  HeatXfer = apex(heatXfer_cospan)

  bespoke_op1_inf_rules = [
  # Rules for avg where tgt is unknown.
  (src_type = :Form0, tgt_type = :infer, replacement_type = :Form1, op = :avg),
  # Rules for avg where src is unknown.
  (src_type = :infer, tgt_type = :Form1, replacement_type = :Form0, op = :avg),
  # Rules for R₀ where tgt is unknown.
  (src_type = :Form0, tgt_type = :infer, replacement_type = :Form0, op = :R₀),
  # Rules for R₀ where src is unknown.
  (src_type = :infer, tgt_type = :Form0, replacement_type = :Form0, op = :R₀),
  # Rules for neg where tgt is unknown.
  (src_type = :Form0, tgt_type = :infer, replacement_type = :Form0, op = :neg),
  (src_type = :Form1, tgt_type = :infer, replacement_type = :Form1, op = :neg),
  (src_type = :Form2, tgt_type = :infer, replacement_type = :Form2, op = :neg),
  # Rules for neg where src is unknown.
  (src_type = :infer, tgt_type = :Form0, replacement_type = :Form0, op = :neg),
  (src_type = :infer, tgt_type = :Form1, replacement_type = :Form1, op = :neg),
  (src_type = :infer, tgt_type = :Form2, replacement_type = :Form2, op = :neg),
  # Rules for Δ where tgt is unknown.
  (src_type = :Form1, tgt_type = :infer, replacement_type = :Form1, op = :Δ),
  # Rules for Δ where src is unknown.
  (src_type = :infer, tgt_type = :Form1, replacement_type = :Form1, op = :Δ),
  # Rules for δ where tgt is unknown.
  (src_type = :Form1, tgt_type = :infer, replacement_type = :Form0, op = :δ),
  # Rules for δ where src is unknown.
  (src_type = :infer, tgt_type = :Form0, replacement_type = :Form1, op = :δ)]

  bespoke_op2_inf_rules = [
  # Rules for / where res is unknown.
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :/),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :infer, replacement_type = :Form1, op = :/),
  (proj1_type = :Form2, proj2_type = :Form2, res_type = :infer, replacement_type = :Form2, op = :/),
  (proj1_type = :Constant, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :/),
  (proj1_type = :Constant, proj2_type = :Form1, res_type = :infer, replacement_type = :Form1, op = :/),
  (proj1_type = :Constant, proj2_type = :Form2, res_type = :infer, replacement_type = :Form2, op = :/),
  (proj1_type = :Form0, proj2_type = :Constant, res_type = :infer, replacement_type = :Form0, op = :/),
  (proj1_type = :Form1, proj2_type = :Constant, res_type = :infer, replacement_type = :Form1, op = :/),
  (proj1_type = :Form2, proj2_type = :Constant, res_type = :infer, replacement_type = :Form2, op = :/),
  # Rules for * where res is unknown.
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :*),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :infer, replacement_type = :Form1, op = :*),
  (proj1_type = :Form2, proj2_type = :Form2, res_type = :infer, replacement_type = :Form2, op = :*),
  (proj1_type = :Constant, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :*),
  (proj1_type = :Constant, proj2_type = :Form1, res_type = :infer, replacement_type = :Form1, op = :*),
  (proj1_type = :Constant, proj2_type = :Form2, res_type = :infer, replacement_type = :Form2, op = :*),
  (proj1_type = :Form0, proj2_type = :Constant, res_type = :infer, replacement_type = :Form0, op = :*),
  (proj1_type = :Form1, proj2_type = :Constant, res_type = :infer, replacement_type = :Form1, op = :*),
  (proj1_type = :Form2, proj2_type = :Constant, res_type = :infer, replacement_type = :Form2, op = :*)]

  infer_types!(HeatXfer, vcat(bespoke_op1_inf_rules, default_op1_type_inference_rules_2D),
    vcat(bespoke_op2_inf_rules, default_op2_type_inference_rules_2D))

  names_types_hx = Set(zip(HeatXfer[:name], HeatXfer[:type]))

  names_types_expected_hx = [
  (Symbol("continuity_advection_•4"), :Form0),
  (:navierstokes_G, :Form1),
  (Symbol("navierstokes_•8"), :Form1),
  (Symbol("navierstokes_•18"), :Form1),
  (Symbol("navierstokes_•3"), :Form1),
  (Symbol("navierstokes_•13"), :Form0),
  (Symbol("navierstokes_•10"), :Form0),
  (Symbol("continuity_diffusion_•2"), :Form1),
  (:continuity_Ṫₐ, :Form0),
  (Symbol("navierstokes_•15"), :Form1),
  (Symbol("continuity_advection_•2"), :Form0),
  (Symbol("navierstokes_•1"), :Form1),
  (Symbol("continuity_diffusion_•1"), :Form1),
  (Symbol("continuity_advection_•3"), :DualForm2),
  (:navierstokes_Ṁ, :Form1),
  (:ρ, :Form0),
  (:T, :Form0),
  (:continuity_diffusion_k, :Constant),
  (Symbol("navierstokes_•7"), :Form1),
  (Symbol("navierstokes_•22"), :DualForm2),
  (Symbol("continuity_diffusion_•3"), :DualForm2),
  (Symbol("navierstokes_•20"), :DualForm2),
  (:continuity_Ṫ, :Form0),
  (:M, :Form1),
  (:navierstokes_three, :Constant),
  (:navierstokes_kᵥ, :Constant),
  (Symbol("navierstokes_•12"), :Form1),
  (:navierstokes_ṗ, :Form0),
  (Symbol("navierstokes_•11"), :Form1),
  (Symbol("navierstokes_•17"), :Form1),
  (Symbol("navierstokes_•6"), :Form1),
  (Symbol("navierstokes_•19"), :Form1),
  (:navierstokes_two, :Constant),
  (:navierstokes_sum_1, :Form1),
  (Symbol("navierstokes_•2"), :Form1),
  (Symbol("navierstokes_•5"), :Form1),
  (Symbol("continuity_advection_•1"), :Form1),
  (Symbol("navierstokes_•4"), :Form1),
  (Symbol("navierstokes_•16"), :Form1),
  (:P, :Form0),
  (:continuity_diffusion_ϕ, :DualForm1),
  (Symbol("navierstokes_•9"), :Form1),
  (:navierstokes_V, :Form1),
  (:continuity_Ṫ₁, :Form0),
  (Symbol("continuity_advection_•5"), :DualForm2),
  (:continuity_advection_V, :Form1),
  (Symbol("navierstokes_•14"), :Form0),
  (Symbol("navierstokes_•21"), :Form0)]

  @test issetequal(names_types_hx, names_types_expected_hx)

  bespoke_op2_res_rules = [
  # Rules for L.
  (proj1_type = :Form1, proj2_type = :DualForm2, res_type = :Form2, resolved_name = :L₀, op = :L),
  (proj1_type = :Form1, proj2_type = :DualForm2, res_type = :DualForm2, resolved_name = :L₀, op = :L),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form1, resolved_name = :L₁′, op = :L)]

  resolve_overloads!(HeatXfer, default_op1_overloading_resolution_rules_2D,
    vcat(bespoke_op2_res_rules, default_op2_overloading_resolution_rules_2D))

  op1s_hx = HeatXfer[:op1]
  op1s_expected_hx = [:d₀, :⋆₁, :dual_d₁, :⋆₀⁻¹, :avg, :R₀, :⋆₀, :⋆₀⁻¹, :neg, :∂ₜ, :avg, :neg, :avg, :Δ₁, :δ₁, :d₀, :d₀, :avg, :d₀, :neg, :avg, :∂ₜ, :⋆₀, :⋆₀⁻¹, :neg, :∂ₜ]
  @test op1s_hx == op1s_expected_hx
  op2s_hx = HeatXfer[:op2]
  op2s_expected_hx = [:*, :/, :/, :L₀, :/, :L₁′, :*, :/, :*, :i₁, :/, :*, :*, :L₀]
  @test op2s_hx == op2s_expected_hx
end
