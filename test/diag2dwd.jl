using Test
using Catlab
using Catlab.Theories
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

@testset "Parsing" begin

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
  js = [Judgement(:C, :Form0, :X), 
        Judgement(:Ċ₁, :Form0, :X),
        Judgement(:Ċ₂, :Form0, :X)
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
  
  sup_js = js = [Judgement(:C, :Form0, :X), 
  Judgement(:ϕ₁, :Form0, :X),
  Judgement(:ϕ₂, :Form0, :X)
  ]
  sup_eqs = [Eq(Var(:ϕ₁), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
         Eq(Var(:ϕ₂), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
         Eq(Tan(Var(:C)), Plus([Var(:ϕ₁), Var(:ϕ₂)]))
  ]
  sup_d = DecaExpr(sup_js, sup_eqs)
  # sup_cset = Decapode(sup_d)
  sup_cset_named = SummationDecapode(sup_d)
  
  
  # Decapodes.compile(diffusion_cset_named, [:C,])
  # Decapodes.compile(test_cset_named, [:C,])
  # Decapodes.compile(sup_cset_named, [:C,])
  
  term(:(∧₀₁(C,V)))

  # No need to parameterize forms over a space (i.e. {X} syntax)
  DiffusionExprBody1 =  quote
    (C, Ċ)::Form0
    ϕ::Form1
    ϕ ==  ∘(k, d₀)(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end
  diffExpr1 = parse_decapode(DiffusionExprBody1)
  ddp1 = SummationDecapode(diffExpr1)
  @test ddp1[:name] == [:C, :Ċ, :ϕ]

  # Support parsing literals.
  DiffusionExprBody2 =  quote
    (C, Ċ)::Form0
    ϕ::Form1
    ϕ ==  2*d₀(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end
  diffExpr2 = parse_decapode(DiffusionExprBody2)
  ddp2 = SummationDecapode(diffExpr2)
  @test ddp2[:name] == [:C, :Ċ, :ϕ, Symbol("•",1), Symbol(2)]

  # Multiply without explicitly giving parentheses.
  DiffusionExprBody3 =  quote
    (C, Ċ)::Form0
    ϕ::Form1
    ϕ ==  4*2*3*d₀(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end
  diffExpr3 = parse_decapode(DiffusionExprBody3)
  ddp3 = SummationDecapode(diffExpr3)
  @test ddp3[:name] == [:C, :Ċ, :ϕ, Symbol("4"), Symbol("2"), Symbol("3"), Symbol("•1"), :mult_1, :mult_2]

  # Variables need not be declared before use.
  DiffusionExprBody4 =  quote
    Ċ::Form0
    ϕ ==  2*d₀(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end
  diffExpr4 = parse_decapode(DiffusionExprBody4)
  ddp4 = SummationDecapode(diffExpr4)
  DiffusionExprBody5 =  quote
    ϕ ==  2*d₀(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end
  diffExpr5 = parse_decapode(DiffusionExprBody5)
  ddp5 = SummationDecapode(diffExpr5)
  @test ddp5[:name] == [:ϕ, :Ċ, Symbol("2"), Symbol("•1"), :C]

  # TVars can be parsed on either side of an equation.
  DiffusionExprBody6 =  quote
    ϕ ==  2*d₀(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    Ċ == ∂ₜ(C)
  end
  diffExpr6 = parse_decapode(DiffusionExprBody6)
  ddp6 = SummationDecapode(diffExpr6)
  @test ddp6[:name] == [:ϕ, :Ċ, Symbol("2"), Symbol("•1"), :C]
  @test ddp6[:incl] == [2]
  DiffusionExprBody7 =  quote
    ϕ ==  2*d₀(C)
    ∂ₜ(C) == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
  end
  diffExpr7 = parse_decapode(DiffusionExprBody7)
  ddp7 = SummationDecapode(diffExpr7)
  @test ddp7[:name] == [:ϕ, :Ċ, Symbol("2"), Symbol("•2"), :C]
  @test ddp7[:incl] == [2]

  # Vars can only be of certain types.
  DiffusionExprBody8 =  quote
    (C)::Foo
    ϕ ==  2*d₀(C)
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    Ċ == ∂ₜ(C)
  end
  diffExpr8 = parse_decapode(DiffusionExprBody8)
  @test_throws ErrorException SummationDecapode(diffExpr8)

  # Multiple equality is not an accepted input
  ParseTest1 = quote
    (A, B, X)::Form0{X}
    A == d(B) == f(X)
  end
  @test_throws ErrorException parse_decapode(ParseTest1)

  # Just noting that the first decapode is denotes a X as an op1
  # while the second is a multiplication between X and F
  ParseTest2_1 = quote
    (A, B, X)::Form0{X}
    A == X(F)
  end
  pt2_1 = SummationDecapode(parse_decapode(ParseTest2_1))
  ParseTest2_2 = quote
    (A, B, X)::Form0{X}
    A == (X)F
  end
  pt2_2 = SummationDecapode(parse_decapode(ParseTest2_2))
  @test pt2_1 != pt2_2

  # Chained Tvars test
  # TODO: Do we want explict support for higher order Tvars?
  ParseTest3 = quote
    D == ∂ₜ(C)
    E == ∂ₜ(D)
  end
  pt3 = SummationDecapode(parse_decapode(ParseTest3))
  @test pt3[:name] == [:D, :E, :C]
  @test pt3[:incl] == [1,2]
  @test pt3[:src] == [3, 1]
  @test pt3[:tgt] == [1, 2]

  dot_rename!(pt3)
  @test pt3[:name] == [Symbol('C'*'\U0307'), Symbol('C'*'\U0307'*'\U0307'), :C]

  # TODO: We should eventually recognize this equivalence
  #= ParseTest4 = quote
    D == D + C
    D + C == C
  end
  pt4 = SummationDecapode(parse_decapode(ParseTest4)) =#

  # Do not rename TVars if they are given a name.
  pt5 = SummationDecapode(parse_decapode(quote
    X::Form0{Point}
    V::Form0{Point}

    k::Constant{Point}

    ∂ₜ(X) == V
    ∂ₜ(V) == -1*k*(X)
  end))

  @test pt5[:name] == [:X, :V, :k, :mult_1, Symbol('V'*'\U0307'), Symbol("-1")]
  dot_rename!(pt5)
  @test pt5[:name] == [:X, Symbol('X'*'\U0307'), :k, :mult_1, Symbol('X'*'\U0307'*'\U0307'), Symbol("-1")]

end
Deca = quote
  (A, B, C)::Form0
end
@testset "Term Construction" begin
    @test term(:(Ċ)) == Var(:Ċ)
    @test_throws ErrorException term(:(∂ₜ{Form0}))
    # @test term(Expr(:ϕ)) == Var(:ϕ)
    @test typeof(term(:(d₀(C)))) == App1
    @test typeof(term(:(∘(k, d₀)(C)))) == AppCirc1
    # @test typeof(term(:(∘(k, d₀)(C,Φ)))) == AppCirc2
    # @test term(:(∘(k, d₀)(C))) == AppCirc1([:k, :d₀], Var(:C)) #(:App1, ((:Circ, :k, :d₀), Var(:C)))
    # @test term(:(∘(k, d₀{X})(C))) == (:App1, ((:Circ, :k, :(d₀{X})), Var(:C)))
    @test_throws MethodError term(:(Ċ == ∘(⋆₀⁻¹{X}, dual_d₁{X}, ⋆₁{X})(ϕ)))
    @test term(:(∂ₜ(C))) == Tan(Var(:C))
    # @test term(:(∂ₜ{Form0}(C))) == App1(:Tan, Var(:C))
end

#= @testset "Recursive Expr" begin
  Recursion = quote
    x::Form0{X}
    y::Form0{X}
    z::Form0{X}

    ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
    y == F(f2(x), ρ(x,z))
  end

  recExpr = parse_decapode(Recursion)
  rdp = SummationDecapode(recExpr)

  @test nparts(rdp, :Var) == 9
  @test nparts(rdp, :TVar) == 1
  @test nparts(rdp, :Op1) == 5
  @test nparts(rdp, :Op2) == 2
  @test nparts(rdp, :Σ) == 1
end =#

@testset "Recursive Expr" begin
  Recursion = quote
    x::Form0{X}
    y::Form0{X}
    z::Form0{X}

    ∂ₜ(z) == f1(x) + ∘(g, h)(y)
    y == F(f2(x), ρ(x,z))
  end

  recExpr = parse_decapode(Recursion)
  rdp = SummationDecapode(recExpr)

  @test nparts(rdp, :Var) == 8
  @test nparts(rdp, :TVar) == 1
  @test nparts(rdp, :Op1) == 4
  @test nparts(rdp, :Op2) == 2
  @test nparts(rdp, :Σ) == 1
end

#= Recursion = quote
  x::Form0{X}
  y::Form0{X}
  z::Form0{X}

  ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
  y == F(f2(x), ρ(x,z))
end =#

Recursion = quote
  x::Form0{X}
  y::Form0{X}
  z::Form0{X}

  ∂ₜ(z) == f1(x) + ∘(g, h)(y)
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
  names_types_expected_1 = Set([(:C, :Form0)])
  @test issetequal(names_types_1, names_types_expected_1)

  # The type of the src of ∂ₜ is inferred.
  Test2 = quote
    C::infer{X}
    ∂ₜ(C) == C
  end
  t2 = SummationDecapode(parse_decapode(Test2))
  t2[only(incident(t2, :C, :name)), :type] = :Form0
  infer_types!(t2)

  names_types_2 = Set(zip(t2[:name], t2[:type]))
  names_types_expected_2 = Set([(:C, :Form0)])
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

  # Special case of a summation with a mix of infer and Literal.
  t16 = @decapode begin
    A == 2 + C + D
  end
  infer_types!(t16)

  names_types_16 = get_name_type_pair(t16)
  names_types_expected_16 = Set([(:A, :Constant), (:C, :Constant), (:D, :Constant), (Symbol("2"), :Literal)])
  @test issetequal(names_types_16, names_types_expected_16)

  # Special case of a summation with a mix of Form, Constant, and infer.
  t17 = @decapode begin
    A::Form0
    C::Constant
    A == C + D
  end
  infer_types!(t17)

  names_types_17 = get_name_type_pair(t17)
  names_types_expected_17 = Set([(:A, :Form0), (:C, :Constant), (:D, :Form0)])
  @test issetequal(names_types_17, names_types_expected_17)

  # Special case of a summation with a mix of infer and Constant.
  t18 = @decapode begin
    C::Constant
    A == C + D
  end
  infer_types!(t18)

  names_types_18 = get_name_type_pair(t18)
  names_types_expected_18 = Set([(:A, :Constant), (:C, :Constant), (:D, :Constant)])
  @test issetequal(names_types_18, names_types_expected_18)
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
    k::Parameter{X}
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
    (two,three,kᵥ)::Parameter{X}
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
  # Rules for avg.
  (src_type = :Form0, tgt_type = :Form1, op_names = [:avg]),
  # Rules for R₀.
  (src_type = :Form0, tgt_type = :Form0, op_names = [:R₀])]

  bespoke_op2_inf_rules = [
  (proj1_type = :Parameter, proj2_type = :Form0, res_type = :Form0, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Parameter, proj2_type = :Form1, res_type = :Form1, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Parameter, proj2_type = :Form2, res_type = :Form2, op_names = [:/, :./, :*, :.*]),

  (proj1_type = :Form0, proj2_type = :Parameter, res_type = :Form0, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Form1, proj2_type = :Parameter, res_type = :Form1, op_names = [:/, :./, :*, :.*]),
  (proj1_type = :Form2, proj2_type = :Parameter, res_type = :Form2, op_names = [:/, :./, :*, :.*])]

  infer_types!(HeatXfer, vcat(bespoke_op1_inf_rules, op1_inf_rules_2D),
    vcat(bespoke_op2_inf_rules, op2_inf_rules_2D))

  names_types_hx = zip(HeatXfer[:name], HeatXfer[:type])

  names_types_expected_hx = [
    (:T, :Form0), (:continuity_Ṫ₁, :Form0), (:continuity_diffusion_ϕ, :DualForm1), (:continuity_diffusion_k, :Parameter), (Symbol("continuity_diffusion_•1"), :DualForm2), (Symbol("continuity_diffusion_•2"), :Form1), (Symbol("continuity_diffusion_•3"), :Form1), (:M, :Form1), (:continuity_advection_V, :Form1), (:ρ, :Form0), (:P, :Form0), (:continuity_Ṫₐ, :Form0), (Symbol("continuity_advection_•1"), :Form0), (Symbol("continuity_advection_•2"), :Form1), (Symbol("continuity_advection_•3"), :DualForm2), (Symbol("continuity_advection_•4"), :Form0), (Symbol("continuity_advection_•5"), :DualForm2), (:continuity_Ṫ, :Form0), (:navierstokes_Ṁ, :Form1), (:navierstokes_G, :Form1), (:navierstokes_V, :Form1), (:navierstokes_ṗ, :Form0), (:navierstokes_two, :Parameter), (:navierstokes_three, :Parameter), (:navierstokes_kᵥ, :Parameter), (Symbol("navierstokes_•1"), :DualForm2), (Symbol("navierstokes_•2"), :Form1), (Symbol("navierstokes_•3"), :Form1), (Symbol("navierstokes_•4"), :Form1), (Symbol("navierstokes_•5"), :Form1), (Symbol("navierstokes_•6"), :Form1), (Symbol("navierstokes_•7"), :Form1), (Symbol("navierstokes_•8"), :Form1), (Symbol("navierstokes_•9"), :Form1), (Symbol("navierstokes_•10"), :Form1), (Symbol("navierstokes_•11"), :Form0), (:navierstokes_sum_1, :Form1), (Symbol("navierstokes_•12"), :Form1), (Symbol("navierstokes_•13"), :Form1), (Symbol("navierstokes_•14"), :Form0), (Symbol("navierstokes_•15"), :Form0), (Symbol("navierstokes_•16"), :Form1), (Symbol("navierstokes_•17"), :Form1), (Symbol("navierstokes_•18"), :Form1), (Symbol("navierstokes_•19"), :Form1), (Symbol("navierstokes_•20"), :Form1), (:navierstokes_sum_2, :Form0), (Symbol("navierstokes_•21"), :DualForm2)]

  @test issetequal(names_types_hx, names_types_expected_hx)

  bespoke_op2_res_rules = [
  # Rules for L.
  (proj1_type = :Form1, proj2_type = :DualForm2, res_type = :Form2, resolved_name = :L₀, op = :L),
  (proj1_type = :Form1, proj2_type = :DualForm2, res_type = :DualForm2, resolved_name = :L₀, op = :L),
  (proj1_type = :Form1, proj2_type = :Form1, res_type = :Form1, resolved_name = :L₁′, op = :L)]

  resolve_overloads!(HeatXfer, op1_res_rules_2D,
    vcat(bespoke_op2_res_rules, op2_res_rules_2D))

  op1s_hx = HeatXfer[:op1]
  op1s_expected_hx = [:d₀, :⋆₁, :dual_d₁, :⋆₀⁻¹, :avg, :R₀, :⋆₀, :⋆₀⁻¹, :neg, :∂ₜ, :avg, :neg, :avg, :Δ₁, :δ₁, :d₀, :d₀, :avg, :d₀, :neg, :avg, :∂ₜ, :⋆₀, :⋆₀⁻¹, :neg, :∂ₜ]
  @test op1s_hx == op1s_expected_hx
  op2s_hx = HeatXfer[:op2]
  op2s_expected_hx = [:*, :/, :/, :L₀, :/, :L₁′, :*, :/, :*, :i₁, :/, :*, :*, :L₀]
  @test op2s_hx == op2s_expected_hx
end

@testset "Compilation Transformation" begin

  # contract_operators does not change the order of op1s.
  Test1 = quote
    (A,B)::Form0{X}
    B == ∘(j,i,h,g,f)(A)
  end
  t1_orig = SummationDecapode(parse_decapode(Test1))
  t1_expanded = expand_operators(t1_orig)
  t1_contracted = contract_operators(t1_expanded)
  @test t1_orig == t1_contracted
  # contract_operators does not mutate its argument.
  @test t1_contracted !== t1_expanded

  # contract_operators works on multiple chains.
  Test2 = quote
    (A,B,C,D)::Form0{X}
    B == ∘(d,d,⋆,d,d)(A)
    D == ∘(d,d,⋆,d,d)(C)
  end
  t2_orig = SummationDecapode(parse_decapode(Test2))
  t2_expanded = expand_operators(t2_orig)
  t2_contracted = contract_operators(t2_expanded)
  #@test t2_orig == t2_contracted
  @test issetequal(
                   [(t2_orig[:src][i], t2_orig[:tgt][i], t2_orig[:op1][i]) for i in parts(t2_orig, :Op1)],
                   [(t2_contracted[:src][i], t2_contracted[:tgt][i], t2_contracted[:op1][i]) for i in parts(t2_contracted, :Op1)])


  # contract_operators "absorbs ∘s".
  t3_orig = SummationDecapode(parse_decapode(quote
    (A,B)::Form0{X}
    B == f(∘(g)(A))
    #D == ∘(g)(f(C)) # Uncomment when parsing supports this
  end))
  t3_expanded = expand_operators(t3_orig)
  t3_contracted = contract_operators(t3_expanded)
  @test t3_contracted == SummationDecapode(parse_decapode(quote
    (A,B)::Form0{X}
    B == ∘(g,f)(A)
  end))

  # contract_operators does not contract through when a var is
  # the target of multiple ops.
  t4_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D)::Form0{X}
    C == f(A)
    C == g(B)
    D == h(C)
  end))
  @test t4_orig == contract_operators(t4_orig)

  # contract_operators does not contract through when a var is a summand.
  t5_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D,E)::Form0{X}
    B == f(A)
    C == g(B)
    E == B + D
  end))
  @test t5_orig == contract_operators(t5_orig)

  # recursive_delete_parents does not change an empty Decapode.
  t6_orig = SummationDecapode(parse_decapode(quote

  end))
  t6_rec_del = recursive_delete_parents(t6_orig, Vector{Int64}())
  @test t6_orig == SummationDecapode{Any, Any, Symbol}()
  # recursive_delete_parents does not mutate its argument.
  @test t6_orig !== t6_rec_del

  # recursive_delete_parents deletes a chain of single-child parents.
  t7_orig = SummationDecapode(parse_decapode(quote
    (A,B,C)::Form0{X}
    B == f(A)
    C == g(B)
  end))
  t7_rec_del = recursive_delete_parents(t7_orig, incident(t7_orig, :C, :name))
  @test t7_rec_del == SummationDecapode{Any, Any, Symbol}()

  # recursive_delete_parents deletes a chain of single-child parents until a
  # multi-child parent is found.
  t8_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D,E)::Form0{X}
    B == f(A)
    C == g(B)
    D == h(C)

    E == i(B)
  end))
  t8_rec_del = recursive_delete_parents(t8_orig, incident(t8_orig, :D, :name))
  @test t8_rec_del == SummationDecapode(parse_decapode(quote
    (A, B, E)::Form0{X}
    B == f(A)

    E == i(B)
  end))

  # recursive_delete_parents deletes a multi-child parent if all its children are
  # deleted.
  t9_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D,E)::Form0{X}
    B == f(A)
    C == g(B)

    D == i(A)
    E == j(D)
  end))
  t9_rec_del = recursive_delete_parents(t9_orig,
                                reduce(vcat, incident(t9_orig, [:C, :E], :name)))
  @test t9_rec_del == SummationDecapode{Any, Any, Symbol}()

  # recursive_delete_parents deletes π₁ and π₂ of an Op2 with missing result.
  t10_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D)::Form0{X}
    C == f(A,B)
    D == g(C)
  end))
  t10_rec_del = recursive_delete_parents(t10_orig, incident(t10_orig, :D, :name))
  @test t10_rec_del == SummationDecapode{Any, Any, Symbol}()

  # recursive_delete_parents only deletes an π₁ or π₂ of an Op2 if it is not used in
  # another Op2.
  t11_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D,E)::Form0{X}
    D == f(A,B)
    E == g(B,C)
  end))
  t11_rec_del = recursive_delete_parents(t11_orig, incident(t11_orig, :D, :name))
  # We should test if Decapodes are isomorphic, but Catlab
  # doesn't do this yet.
  #@test t11_rec_del == SummationDecapode(parse_decapode(quote
  #  (B,C,E)::Form0{X}
  #  E == g(B,C)
  #end))
  @test t11_rec_del == @acset SummationDecapode{Any, Any, Symbol} begin
    Var  = 3
    Op2  = 1

    type = [:Form0, :Form0, :Form0]
    name = [:E, :B, :C]

    proj1 = [2]
    proj2 = [3]
    res   = [1]
    op2   = [:g]
  end

  # If π₁ and π₂ of an Op2 point to the same Var, that Var is used in some
  # other computation, but this Op2 has deleted result, then delete the Op2,
  # but not the Var.
  t12_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D)::Form0{X}
    C == f(A,A)
    D == g(A,B)
  end))
  t12_rec_del = recursive_delete_parents(t12_orig, incident(t12_orig, :C, :name))
  @test t12_rec_del == SummationDecapode(parse_decapode(quote
    (A,B,D)::Form0{X}
    D == g(A,B)
  end))

  # recursive_delete_parents does not delete a Var if it is used in a summation.
  t13_orig = SummationDecapode(parse_decapode(quote
    (A,B,C,D,E)::Form0{X}
    D == A + B + C

    E == f(A)
  end))
  t13_rec_del = recursive_delete_parents(t13_orig, incident(t13_orig, :E, :name))
  @test t13_rec_del == SummationDecapode(parse_decapode(quote
    (A,B,C,D)::Form0{X}
    D == A + B + C
  end))

  # Test recursive_delete_parents on a real physics.
  Veronis = SummationDecapode(parse_decapode(quote
    B::DualForm0{X}
    E::Form1{X}
    σ::Form1{X}
    J::Form1{X}
    (negone,c,ε₀)::Constant{X}

    ∂ₜ(E) == ((negone * (J - σ .* E))./ε₀) + ((c * c).*(⋆₁⁻¹(d₀(B))))

    ∂ₜ(B) == ⋆₂(d₁(E))
  end))
  Veronis_rec_del = recursive_delete_parents(Veronis, incident(Veronis, :Ḃ, :name))
  @test Veronis_rec_del == SummationDecapode(parse_decapode(quote
    B::DualForm0{X}
    E::Form1{X}
    σ::Form1{X}
    J::Form1{X}
    (negone,c,ε₀)::Constant{X}

    ∂ₜ(E) == ((negone * (J - σ .* E))./ε₀) + ((c * c).*(⋆₁⁻¹(d₀(B))))
  end))

  # Test recursive_delete_parents on a real physics.
  Veronis = SummationDecapode(parse_decapode(quote
    B::DualForm0{X}
    E::Form1{X}
    σ::Form1{X}
    J::Form1{X}
    (negone,c,ε₀)::Constant{X}

    ∂ₜ(E) == ((negone * (J - σ .* E))./ε₀) + ((c * c).*(⋆₁⁻¹(d₀(B))))

    ∂ₜ(B) == ⋆₂(d₁(E))
  end))
  Veronis_rec_del = recursive_delete_parents(Veronis, incident(Veronis, :Ė, :name))
  # We should test if Decapodes are isomorphic, but Catlab
  # doesn't do this yet.
  #@test Veronis_rec_del == SummationDecapode(parse_decapode(quote
  #  B::DualForm0{X}
  #  E::Form1{X}

  #  ∂ₜ(B) == ⋆₂(d₁(E))
  #end))
  @test Veronis_rec_del ==  @acset SummationDecapode{Any, Any, Symbol} begin
    Var  = 4
    TVar = 1
    Op1  = 3

    type = [:DualForm0, :Form1, :infer, :infer]
    name = [:B, :E, :sum_1, :Ḃ]
    
    incl = [4]

    src = [3, 1, 2]
    tgt = [4, 4, 3]

    op1 = [:⋆₂, :∂ₜ, :d₁]
  end
end

@testset "ASCII & Vector Calculus Operators" begin
  # Test ASCII to Unicode conversion on an Op2.
  t1 = @decapode begin
    A == wedge(C, D)
  end
  unicode!(t1)

  op2s_1 = Set(t1[:op2])
  op2s_expected_1 = Set([:∧])
  @test issetequal(op2s_1, op2s_expected_1)

  # Test ASCII to Unicode conversion with multiple occurences of the same operator.
  t2 = @decapode begin
    A == wedge(C, D)
    B == wedge(E, F)
  end
  unicode!(t2)

  op2s_2 = Set(t2[:op2])
  op2s_expected_2 = Set([:∧])
  @test issetequal(op2s_2, op2s_expected_2)

  # Test ASCII to Unicode conversion works with composed operators after expansion.
  t3 = @decapode begin
    A == ∘(star, lapl, star_inv)(B)
  end
  t3 = expand_operators(t3)
  unicode!(t3)

  op1s_3 = Set(t3[:op1])
  op1s_expected_3 = Set([:⋆,:Δ,:⋆⁻¹])
  @test issetequal(op1s_3, op1s_expected_3)

  # Test ASCII tangent operator identifies a TVar.
  t4 = @decapode begin
    A == dt(B)
  end
  @test nparts(t4, :TVar) == 1

  # Test vec_to_dec! on a single operator.
  t5 = @decapode begin
    A == div(B)
  end
  vec_to_dec!(t5)

  op1s_5 = Set(t5[:op1])
  op1s_expected_5 = Set([[:⋆,:d,:⋆]])
  @test issetequal(op1s_5, op1s_expected_5)

  # Test divergence of gradient is the Laplacian.
  t6 = @decapode begin
    A == ∘(grad, div)(B)
  end
  t6 = expand_operators(t6)
  vec_to_dec!(t6)
  t6 = contract_operators(t6)
  @test only(t6[:op1]) == [:d,:⋆,:d,:⋆]

  # Test curl of curl is a vector Laplacian.
  t7 = @decapode begin
    A == ∘(∇x, ∇x)(B)
  end
  t7 = expand_operators(t7)
  vec_to_dec!(t7)
  t7 = contract_operators(t7)
  @test only(t7[:op1]) == [:d,:⋆,:d,:⋆]

  # Test advection is divergence of wedge product.
  t8 = @decapode begin
    A == adv(B, C)
  end
  vec_to_dec!(t8)

  @test only(t8[:op1]) == [:⋆,:d,:⋆]
  @test only(t8[:op2]) == :∧

  # Test multiple advections.
  t9 = @decapode begin
    A == adv(B, C)
    D == adv(E, F)
  end
  vec_to_dec!(t9)

  t9_expected = @decapode begin
    A == ∘(⋆,d,⋆)(B ∧ C)
    D == ∘(⋆,d,⋆)(E ∧ F)
  end
  t9_expected[incident(t9_expected, Symbol("•1"), :name), :name] = Symbol("•_adv_1")
  t9_expected[incident(t9_expected, Symbol("•2"), :name), :name] = Symbol("•_adv_2")
  @test is_isomorphic(t9, t9_expected)

end

