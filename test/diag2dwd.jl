using Test
using Catlab
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Present
using Catlab.Programs
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
        V::Form1{X}
        ϕ::Form1{X}

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
    advdiffdp = SummationDecapode(advdiff)
    @test nparts(advdiffdp, :Var) == 6
    @test nparts(advdiffdp, :TVar) == 1
    @test nparts(advdiffdp, :Op1) == 4
    @test nparts(advdiffdp, :Op2) == 1
    @test nparts(advdiffdp, :Σ) == 1
    @test nparts(advdiffdp, :Summand) == 2
end

@testset "Decapode Composition" begin

  """      function unique_by!(acset, column_names::Vector{Symbol})

  Given column names from the same table, remove duplicate rows.

  WARNING: This function does not check if other tables index into the one
  given. Removal of rows is performed with prejudice.

  See also: [`unique_by`](@ref).

  # Examples
  ```julia-repl
  julia> unique_by!(parallel_arrows(Graph, 123), :E, [:src,:tgt]) == parallel_arrows(Graph, 1)
  true
  ```
  """
  function unique_by!(acset, table::Symbol, columns::Vector{Symbol})
    # TODO: Declarative CT methods are prefered to imperative index arithmetic.
    # Replace this.
    rows = mapreduce(x -> acset[x], hcat, columns)
    rem_parts!(acset, table,
      setdiff(parts(acset, table),
        unique(i -> rows[i,:], eachindex(eachrow(rows)))))
    return acset
  end

  """      function unique_by(acset, column_names::Vector{Symbol})

  Given column names from the same table, return a copy of the acset with
  duplicate rows removed. Removal of rows is performed with prejudice.

  WARNING: This function does not check if other tables index into the one
  given. Removal of rows is performed with prejudice.

  See also: [`unique_by!`](@ref).

  # Examples
  ```julia-repl
  julia> unique_by(parallel_arrows(Graph, 123), :E, [:src,:tgt]) == parallel_arrows(Graph, 1)
  true
  ```
  """
  function unique_by(acset, table::Symbol, columns::Vector{Symbol})
    acset_copy = copy(acset)
    unique_by!(acset_copy, table, columns)
  end

  """    function type_check_decapodes_composition(decapodes_vars, relation, local_ports)

  Check that the types of all Vars connected by the same junction match.

  This function only throws an error on the first type mismatch found.
  """
  function type_check_decapodes_composition(decapodes_vars, relation, local_ports)
    r = relation
    decs = first.(decapodes_vars)
    vars = last.(decapodes_vars)
    for j ∈ junctions(r)
      # Get the type of the first variable attached to this junction.
      P = ports_with_junction(r, j)
      p₁ = first(P)
      b₁ = r[p₁, :box]
      lp₁ = local_ports[p₁]
      type₁ = decs[b₁][lp₁, :type]
      # Compare the first type to the rest of the types.
      for p ∈ rest(P, 2)
        b = r[p, :box]
        lp = local_ports[p]
        symbol_name = vars[b][lp]
        var = only(incident(decs[b], symbol_name, :name))
        type = decs[b][var, :type]
        # TODO: We use == here because we assume a type is a Symbol, like
        # :Form0, :Form1. Will a type ever be something we should compare using
        # isequal instead?
        type == type₁ || let
          var_name =  decs[b ][var, :name]
          var_name₁ = decs[b₁][var, :name]
          decapode_name =  r[b,  :name]
          decapode_name₁ = r[b₁, :name]
          error("The type of $(var_name), $(type), in decapode "*
            "\"$(decapode_name)\" does not match the type of $(var_name₁), "*
            "$(type₁), in decapode \"$(decapode_name₁)\". "*
            "(Also, check that the order of the decapodes you supplied "*
            "matches the the order you specified in the relation.)")
        end
      end
    end
  end

  OpenNamedDecapodeOb, OpenNamedDecapode = OpenACSetTypes(NamedDecapode, :Var)

  # TODO: This does not work:
  # function OpenNamedDecapode(relation, decapode, box)
  function MakeOpenNamedDecapode(relation::RelationDiagram,
    decapode::NamedDecapode{Any, Any, Symbol}, box)
    P = ports(relation, box)
    J = relation[P, :junction]
    V = relation[J, :variable]
    FinFunctions = map(V) do v
      FinFunction(incident(decapode, v, :name), nparts(decapode, :Var))
    end
    OpenNamedDecapode{Any, Any, Symbol}(decapode, FinFunctions...)
  end

  """    function compose_decapodes(decapodes_vars::Vector{Tuple{NamedDecapode{Any, Any, Symbol}, Vector{Symbol}}}, relation::RelationDiagram)

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
  end;

  julia> compose_decapodes([(Diffusion, [:C, :ϕ]), (Advection, [:C, :ϕ, :V]),
    (Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])], compose_diff_adv);
  ```
  """
  function compose_decapodes(decapodes_vars::Vector{Tuple{NamedDecapode{
    Any, Any, Symbol}, Vector{Symbol}}}, relation::RelationDiagram)
    # TODO: Calling this function "compose" is something of a misnomer.  e.g.
    # You could very well give a vector containing a single decapode. Then the
    # effect of applying the relation would be something like renaming
    # variables. So this function should really be called oapply_decapodes (or
    # just oapply if you are willing to import oapply from
    # Catlab.WiringDiagrams).
    r = relation
    copies = @. copy(first(decapodes_vars))

    # Check that the number of decapodes given matches the number of boxes in
    # the relation.
    num_boxes = nboxes(r)
    num_decapodes = length(decapodes_vars)
    # TODO: Should this be an ArgumentError?
    num_boxes == num_decapodes || error(
      "$(num_boxes) decapodes were specified in the relation but only "*
      "$(num_decapodes) were given.")

    # Check that the number of variables given in the relation is the same as
    # the number of symbols in the corresponding vector of Vars.
    # TODO: Should this be an ArgumentError?
    for b ∈ boxes(r)
      # Note: This only returns the first length mismatch found.
      num_junctions = length(incident(r, b, :box))
      num_symbols = length(decapodes_vars[b][2])
      num_junctions == num_symbols || let decapode_name = r[b,  :name]
        error("Decapode \"$(decapode_name)\" needs $(num_junctions) "*
        "variables, but $(num_symbols) were specified.")
      end
    end

    ## Determine the mapping of global ports to local ports.
    # In a RelationDiagram, this is baked into the order of rows in the Port
    # table.
    # This is a column that one could hcat to the Ports table.
    local_ports = [lp for b=boxes(r) for lp=eachindex(ports(r, b))]

    ## Check that types of variables connected by the same junction match.
    # "Do typechecking."
    type_check_decapodes_composition(decapodes_vars, relation, local_ports)

    ## Do namespacing.
    # Append each Var name with the name @relation gave the decapode.
    for b ∈ boxes(r)
      box_name = r[b, :name]
      for v ∈ parts(copies[b], :Var)
        var_name = copies[b][v, :name]
        copies[b][v, :name] = Symbol(box_name, '_', var_name)
      end
    end

    ## Write over the name fields to be what was specified by @relation. (oapply
    # cannot combine objects whose attributes are not equal.)
    for p ∈ ports(r)
      b = r[p, :box]
      j = r[p, :junction]
      lp = local_ports[p]
      # Get the index of the row with this name in the Var.
      symbol_name = decapodes_vars[b][2][lp]
      name = Symbol((r[:name][b]), '_', symbol_name)
      # Note: only is not necessary but is a useful check the decapode is
      # well-formed. If we ever want e.g. X:Form0 and X:Form1 in a single
      # decapode, this will need refactoring.
      var = only(incident(copies[b], name, :name))
      copies[b][var, :name] = r[j, :variable]
    end

    ## Compose
    apex(oapply(relation, map(boxes(r)) do b
      MakeOpenNamedDecapode(r, copies[b], b)
    end))
  end

  function compose_decapodes(decapode_vars, relation::RelationDiagram)
    compose_decapodes([decapode_vars], relation)
  end

  # This is the example from the "Overview" page in the docs.
  DiffusionExprBody =  quote
    C::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ ==  ∘(k, d₀)(C)
  end
  AdvectionExprBody = quote
    C::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}

    ϕ == ∧₀₁(C,V)
  end
  SuperpositionExprBody = quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}

    ϕ == ϕ₁ + ϕ₂
    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
    ∂ₜ(C) == Ċ
  end

  difExpr = parse_decapode(DiffusionExprBody)
  Diffusion = NamedDecapode(difExpr)

  advExpr = parse_decapode(AdvectionExprBody)
  Advection = NamedDecapode(advExpr)

  supExpr = parse_decapode(SuperpositionExprBody)
  Superposition = NamedDecapode(supExpr)

  compose_diff_adv = @relation (C,V) begin
    diffusion(C, ϕ₁)
    advection(C, ϕ₂, V)
    superposition(ϕ₁, ϕ₂, ϕ, C)
  end

  decapodes_vars = [
    (Diffusion, [:C, :ϕ]),
    (Advection, [:C, :ϕ, :V]),
    (Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])]

  dif_adv_sup = compose_decapodes(decapodes_vars, compose_diff_adv)
  dif_adv_sup_expected = @acset NamedDecapode{Any, Any, Symbol} begin
    Var = 6
    type = [:Form0, :Form1, :Form1, :Form1, :infer, :Form1]
    name = [:C, :ϕ₁, :V, :ϕ₂, :superposition_Ċ, :ϕ]

    TVar = 1
    incl = [5]

    Op1 = 3
    src = [1,6,1]
    tgt = [2,5,5]
    op1 = [[:k, :d₀], [:⋆₀⁻¹, :dual_d₁, :⋆₁], :∂ₜ]

    Op2 = 2
    proj1 = [1,2]
    proj2 = [3,4]
    res = [4,6]
    op2 = [:∧₀₁, :+]
  end
  @test dif_adv_sup == dif_adv_sup_expected

  # Test some other permutation of the symbols yields the same decapode.
  compose_diff_adv = @relation (C,V) begin
    diffusion(C, ϕ₁)
    advection(V, ϕ₂, C)
    superposition(ϕ₁, ϕ₂, C, ϕ)
  end
  decapodes_vars = [
    (Diffusion, [:C, :ϕ]),
    (Advection, [:V, :ϕ, :C]),
    (Superposition, [:ϕ₁, :ϕ₂, :C, :ϕ])]
  dif_adv_sup = compose_decapodes(decapodes_vars, compose_diff_adv)
  @test dif_adv_sup == dif_adv_sup_expected

  # Test that Op2s are properly de-duplicated.
  AdvectionExprBody = quote
    C::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ == ∧₀₁(C,V)
  end
  advExpr = parse_decapode(AdvectionExprBody)
  Advection = NamedDecapode(advExpr)
  self_adv = @relation () begin
    advection(C,V,ϕ)
    advection(C,V,ϕ)
  end
  adv_adv = [
   (Advection, [:C,:V,:ϕ]),
   (Advection, [:C,:V,:ϕ])]
  adv_adv_comp = compose_decapodes(adv_adv, self_adv)
  # De-duplicate Op1s.
  unique_by!(adv_adv_comp, :Op1, [:src, :tgt, :op1])
  # De-duplicate Op2s.
  unique_by!(adv_adv_comp, :Op2, [:proj1, :proj2, :res, :op2])
  adv_adv_comp_expected = @acset NamedDecapode{Any, Any, Symbol} begin
    Var = 3
    type = [:Form0, :Form1, :Form1]
    name = [:C, :V, :ϕ]
    Op2 = 1
    proj1 = [1]
    proj2 = [2]
    res = [3]
    op2 = [:∧₀₁]
  end
  @test adv_adv_comp == adv_adv_comp_expected

  # Test that a relation of only one decapode is handled properly.
  AdvectionExprBody =  quote
    C::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ == ∧₀₁(C,V)
  end
  advExpr = parse_decapode(AdvectionExprBody)
  Advection = NamedDecapode(advExpr)
  adv_relation = @relation () begin
    advection(C,V,ϕ)
  end
  adv_comp = compose_decapodes((Advection, [:C,:V,:ϕ]), adv_relation)
  adv_comp_expected = @acset NamedDecapode{Any, Any, Symbol} begin
    Var = 3
    type = [:Form0, :Form1, :Form1]
    name = [:C, :V, :ϕ]
    Op2 = 1
    proj1 = [1]
    proj2 = [2]
    res = [3]
    op2 = [:∧₀₁]
  end
  @test adv_comp == adv_comp_expected

  # Simplest possible decapode relation.
  TrivialExprBody = quote
    H::Form0{X}
  end
  trivalExpr = parse_decapode(TrivialExprBody)
  Trivial = NamedDecapode(trivalExpr)
  trivial_relation = @relation () begin
    trivial(H)
  end
  trivial_comp = compose_decapodes((Trivial, [:H]), trivial_relation)
  @test trivial_comp == Trivial

end

advdiff = parse_decapode(AdvDiff)
advdiffdp = SummationDecapode(advdiff)
to_graphviz(advdiffdp)

#AdvDiff = quote
#    C::Form0{X}
#    Ċ::Form0{X}
#    V::Form1{X}
#    ϕ::Form1{X}
#    ϕ₁::Form1{X}
#    ϕ₂::Form1{X}
#
#    # Fick's first law
#    ϕ₁ ==  (k ∘ d₀)(C)
#    ϕ₂ == ∧₀₁(C,V)
#    ϕ == ϕ₁ + ϕ₂
#    # Diffusion equation
#    ∂ₜ(C) == Ċ
#    Ċ == ∘(⋆₀⁻¹, dual_d₁, ⋆₁)(ϕ)
#end
#
#advdiff = parse_decapode(AdvDiff)
#advdiffdp = NamedDecapode(advdiff)
#to_graphviz(advdiffdp)
#
#compile(advdiffdp, [:C, :V])
