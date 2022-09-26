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
       Eq(Tan(Var(:C)), Plus(Var(:Ċ₁), Var(:Ċ₂)))
]
diffusion_d = DecaExpr(js, eqs)
diffusion_cset = Decapode(diffusion_d)
diffusion_cset_named = NamedDecapode(diffusion_d)
# A test with expressions on LHS (i.e. temporary variables must be made)
# TODO: we should be intelligent and realize that the LHS of the first two
# equations are the same and so can share a new variable
eqs = [Eq(Plus(Var(:Ċ₁), Var(:Ċ₂)), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Plus(Var(:Ċ₁), Var(:Ċ₂)), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus(Var(:Ċ₁), Var(:Ċ₂)))    
]
test_d = DecaExpr(js, eqs)
test_cset = Decapode(test_d)
test_cset_named = NamedDecapode(test_d)

# TODO: Write tests for recursive expressions

all(isassigned(test_cset_named[:name], i) for i in parts(test_cset_named,:Var))

sup_js = js = [Judge(Var(:C), :Form0, :X), 
Judge(Var(:ϕ₁), :Form0, :X),
Judge(Var(:ϕ₂), :Form0, :X)
]
sup_eqs = [Eq(Var(:ϕ₁), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :k, :d₀], Var(:C))),
       Eq(Var(:ϕ₂), AppCirc1([:⋆₀⁻¹, :dual_d₁, :⋆₁, :d₀], Var(:C))),
       Eq(Tan(Var(:C)), Plus(Var(:ϕ₁), Var(:ϕ₂)))    
]
sup_d = DecaExpr(sup_js, sup_eqs)
sup_cset = Decapode(sup_d)
sup_cset_named = NamedDecapode(sup_d)


compile(diffusion_cset_named, [:C,])
compile(test_cset_named, [:C,])
compile(sup_cset_named, [:C,])

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
  rdp = NamedDecapode(recExpr)
  show(rdp)

  @test nparts(rdp, :Var) == 9
  @test nparts(rdp, :TVar) == 1
  @test nparts(rdp, :Op1) == 5
  @test nparts(rdp, :Op2) == 3
end
Recursion = quote
  x::Form0{X}
  y::Form0{X}
  z::Form0{X}

  ∂ₜ(k(z)) == f1(x) + ∘(g, h)(y)
  y == F(f2(x), ρ(x,z))
end

recExpr = parse_decapode(Recursion)
rdp = NamedDecapode(recExpr)

@testset "Diffusion Diagram" begin
    DiffusionExprBody =  quote
        C::Form0{X}
        Ċ::Form0{X}
        ϕ::Form1{X}
    
        # Fick's first law
        ϕ ==  ∘(k, d₀)(C)
        # Diffusion equation
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    diffExpr = parse_decapode(DiffusionExprBody)
    ddp = NamedDecapode(diffExpr)
    to_graphviz(ddp)

    @test nparts(ddp, :Var) == 3
    @test nparts(ddp, :TVar) == 1
    @test nparts(ddp, :Op1) == 3
    @test nparts(ddp, :Op2) == 0
end


@testset "Advection Diagram" begin
    Advection = quote
        C::Form0{X}
        V::Form1{X}
        ϕ::Form1{X}

        ϕ == ∧₀₁(C,V)
    end

    advdecexpr = parse_decapode(Advection)
    advdp = NamedDecapode(advdecexpr)
    @test nparts(advdp, :Var) == 3
    @test nparts(advdp, :TVar) == 0
    @test nparts(advdp, :Op1) == 0
    @test nparts(advdp, :Op2) == 1
end

@testset "Superposition Diagram" begin
    Superposition = quote
        C::Form0{X}
        Ċ::Form0{X}
        ϕ::Form1{X}
        ϕ₁::Form1{X}
        ϕ₂::Form1{X}

        ϕ == ϕ₁ + ϕ₂
        Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
        ∂ₜ(C) == Ċ
    end

    superexp = parse_decapode(Superposition)
    supdp = NamedDecapode(superexp)
    @test nparts(supdp, :Var) == 5
    @test nparts(supdp, :TVar) == 1
    @test nparts(supdp, :Op1) == 2
    @test nparts(supdp, :Op2) == 1
end

@testset "AdvectionDiffusion Diagram" begin
    AdvDiff = quote
        C::Form0{X}
        Ċ::Form0{X}
        V::Form1{X}
        ϕ::Form1{X}
        ϕ₁::Form1{X}
        ϕ₂::Form1{X}
    
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
    advdiffdp = NamedDecapode(advdiff)
    @test nparts(advdiffdp, :Var) == 6
    @test nparts(advdiffdp, :TVar) == 1
    @test nparts(advdiffdp, :Op1) == 4
    @test nparts(advdiffdp, :Op2) == 2
end

@testset "Compose via Structured Cospans" begin
  """
      function compose_decapodes(decapodes_vars, relation::RelationDiagram)

  Compose a list of decapodes as specified by the given relation diagram.

  The decapodes must be given in the same order as they were specified in the
  relation.

  State variables (such as the (C,V) given in the head of the following
  @relation) do not affect the result of a composition.

  # Examples
  ```julia-repl
  julia> compose_diff_adv = @relation (C,V) begin
    diffusion(C, ϕ₁)
    advection(C, ϕ₂, V)
    superposition(ϕ₁, ϕ₂, ϕ, C)
  end
  julia> compose_decapodes([(Diffusion, [:C, :ϕ]), (Advection, [:C, :ϕ, :V]),
    (Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])], compose_diff_adv)
  ````
  """
  function compose_decapodes(decapodes_vars, relation::RelationDiagram)
    r = relation
    copies = map(d -> copy(d[1]), decapodes_vars)
    # TODO: We should also check that the number of variables given in the
    # relation is the same as the number of Vars in the corresponding
    # decapodes.
    # Step -3: Determine the mapping of global ports to local ports. (i.e. for
    # each (box,junction) pair, determine whether this corresponds to the 1st
    # Var of that box, or the 2nd Var of that box, etc.)
    # tuples (box, indices of rows with that box):
    box_rows = [(b, findall(==(b), r[:box])) for b in unique(r[:box])]
    # Decrement each index by the offset of the first index of the group from
    # row 1.
    local_ports = [map(y -> y-idxs[begin]+1, idxs) for (b, idxs) in box_rows]
    # Put into a column that one could append to the Ports table.
    local_ports = local_ports |> flatten |> collect

    # Step -2: Check that the types of all Vars connected by the same
    # junction are the same.
    foreach(unique(r[:junction])) do j
      j_idxs = findall(==(j), r[:junction])
      first_var_idx = copies[r[:box][j_idxs]][begin]
      first_lp_idx = local_ports[j_idxs][begin]
      first_type = first_var_idx[:type][first_lp_idx]
      foreach(r[:box][j_idxs], local_ports[j_idxs]) do b_idx, lp_idx
        # TODO: The `infer` type will never be here once the decapodes parsing
        # is finished being refactored. So we don't check for it here.
        # TODO: Make this error message more informative via string interpolation.
        copies[b_idx][:type][lp_idx] == first_type || error("Types do not match.")
      end
    end

    # Step -1: Make copies of the decapodes, and write over the name fields
    # to be what was specified by @relation.
    foreach(r[:box], r[:junction], local_ports) do b_idx, j_idx, lp_idx
      copies[b_idx][:name][lp_idx] = relation[:variable][j_idx]
    end

    # Step 1: Convert from name to index in the Var tables. (TODO: How much
    # of this bookkeeping does @relation already handle?)
    # TODO: We should create the FinFunctions inside the composition loop?
    # TODO: The order in which you compose cospans matters. Take care when
    # creating FinFunctions.

    # Step 2: Start composing.
    NamedDecapodeOVOb, NamedDecapodeOV = OpenACSetTypes(NamedDecapode, :Var)
    # TODO: dec is a placeholder variable so we can do steps 3 and 4.
    dec = nothing

    # TODO: Perhaps split the de-duplication steps into a single function,
    # or create a function that generalizes de-duplication over all ACSets
    # and all Objects.
    # Step 3: De-duplicate Op1s.
    # In SQL: Select DISTINCT src, tgt, op1 FROM Op1;
    op1_tuples = zip(dec[:src], dec[:tgt], dec[:op1]) |> collect
    rem_parts!(dec, :Op1,
      setdiff(parts(dec, :Op1),
        unique(i -> op1_tuples[i], 1:length(op1_tuples))))

    # Step 4: De-duplicate Op2s.
    # In SQL: Select DISTINCT proj1, proj2, res, op2 FROM Op2;
    op2_tuples = zip(dec[:proj1], dec[:proj2], dec[:res], dec[:op2]) |> collect
    rem_parts!(dec, :Op2,
      setdiff(parts(dec, :Op2),
        unique(i -> op2_tuples[i], 1:length(op2_tuples))))

    # TODO: In case we have to de-duplicate tvars, this is the code. The zip |>
    # collect idiom is just use to mirror the prior calls de-dups.
    ## Step 5: De-duplicate TVars.
    ## In SQL: Select DISTINCT incl FROM TVar;
    #tvar_tuples = zip(dec[:incl]) |> collect
    #rem_parts!(dec, :TVar,
    #  setdiff(parts(dec, :TVar),
    #    unique(i -> tvar_tuples[i], 1:length(tvar_tuples))))
  end

  DiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}
  
    # Fick's first law
    ϕ ==  ∘(k, d₀)(C)
    # Diffusion equation
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end
  AdvectionExprBody =  quote
    C::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}

    ϕ == ∧₀₁(C,V)
  end

  diffExpr = parse_decapode(DiffusionExprBody)
  ddp = NamedDecapode(diffExpr)

  advExpr = parse_decapode(AdvectionExprBody)
  adp = NamedDecapode(advExpr)

  OpenNamedDecapodeOb, OpenNamedDecapode = OpenACSetTypes(NamedDecapode, :Var)

  # TODO: Rewrite this test by calling the compose_decapodes function when finished.
  ddpov = OpenNamedDecapode{Any, Any, Symbol}(ddp, FinFunction([1,3], 3), FinFunction([1,3], 3))
  adpov = OpenNamedDecapode{Any, Any, Symbol}(adp, FinFunction([1,3], 3), FinFunction([1,3], 3))
  advdiff = apex(compose(ddpov, adpov))

  advdiff_expected = @acset NamedDecapode{Any, Any, Symbol} begin
    Var = 4
    type = [:Form0, :infer, :Form1, :Form1]
    name = [:C, :Ċ, :ϕ, :V]
    
    TVar = 1
    incl = [2]
    
    Op1 = 3
    src = [1,3,1]
    tgt = [3,2,2]
    op1 = [[:k, :d₀], [:⋆₀⁻¹, :dual_d₁, :⋆₁], :∂ₜ]
  
    Op2 = 1
    proj1 = [1]
    proj2 = [4]
    res = [3]
    op2 = [:∧₀₁]
  end
  # TODO: This should really be a test for equality between canonical forms
  # (since order within tables does not matter save for preserving mappings
  # between tables.)
  @test advdiff == advdiff_expected

  # TODO: Add a test that checks that there are no duplicate records in the Op1 table.

  # Test that there are no duplicates in the Op2 table.
  # TODO: Rewrite this test by calling the compose_decapodes function when finished.
  adpov_self = OpenNamedDecapode{Any, Any, Symbol}(adp, FinFunction([1,2,3], 3), FinFunction([1,2,3], 3))
  adpov_self_compd = apex(compose(adpov_self, adpov_self))
  op2_tuples = zip(adpov_self_compd[:proj1], adpov_self_compd[:proj2], adpov_self_compd[:res], adpov_self_compd[:op2]) |> collect
  rem_parts!(adpov_self_compd, :Op2,
      setdiff(parts(adpov_self_compd, :Op2),
          unique(i -> op2_tuples[i], 1:length(op2_tuples))))
  #┌─────┬───────┬──────┐
  #│ Var │  type │ name │
  #├─────┼───────┼──────┤
  #│   1 │ Form0 │    C │
  #│   2 │ Form1 │    V │
  #│   3 │ Form1 │    ϕ │
  #└─────┴───────┴──────┘
  #┌─────┬───────┬───────┬─────┬─────┐
  #│ Op2 │ proj1 │ proj2 │ res │ op2 │
  #├─────┼───────┼───────┼─────┼─────┤
  #│   2 │     1 │     2 │   3 │ ∧₀₁ │
  #└─────┴───────┴───────┴─────┴─────┘
  adpov_self_expected = @acset NamedDecapode{Any, Any, Symbol} begin
    Var = 3
    type = [:Form0, :Form1, :Form1]
    name = [:C, :V, :ϕ]
  
    Op2 = 1
    proj1 = [1]
    proj2 = [2]
    res = [3]
    op2 = [:∧₀₁]
  end
  # TODO: This should really be a test for equality between canonical forms
  # (since order within tables does not matter save for preserving mappings
  # between tables.)
  @test adp == adpov_self_expected
end


AdvDiff = quote
    C::Form0{X}
    Ċ::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}

    # Fick's first law
    ϕ₁ ==  (k ∘ d₀)(C)
    ϕ₂ == ∧₀₁(C,V)
    ϕ == ϕ₁ + ϕ₂
    # Diffusion equation
    ∂ₜ(C) == Ċ
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
end

advdiff = parse_decapode(AdvDiff)
advdiffdp = NamedDecapode(advdiff)
to_graphviz(advdiffdp)

compile(advdiffdp, [:C, :V])
