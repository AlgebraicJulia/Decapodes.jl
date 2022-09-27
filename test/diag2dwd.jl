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

@testset "Decapode Composition" begin

  """
      function unique_by!(acset, column_names::Vector{Symbol})

  Given column names from the same table, remove duplicate rows.
  
  WARNING: This function does not check if other tables index into the one
  given. Removal of rows is performed with prejudice.

  See also: [`unique_by`](@ref).

  # Examples
  ```julia-repl
  julia> g = @acset Graph begin
    V = 2
    E = 3
    src = [1,1,2]
    tgt = [1,1,2]
  end
  julia> unique_by!(g, :E, [:src, :tgt])
  ``` 
  """
  function unique_by!(acset, table::Symbol, columns::Vector{Symbol})
    rows = columns .|> (x -> acset[x]) |> (x -> zip(x...)) |> collect
    rem_parts!(acset, table,
      setdiff(parts(acset, table),
        unique(i -> rows[i], eachindex(rows))))
    return acset
  end 

  """
      function unique_by(acset, column_names::Vector{Symbol})

  Given column names from the same table, return a copy of the acset with
  duplicate rows removed. Removal of rows is performed with prejudice.

  WARNING: This function does not check if other tables index into the one
  given. Removal of rows is performed with prejudice.

  See also: [`unique_by!`](@ref).

  # Examples
  ```julia-repl
  julia> g = @acset Graph begin
    V = 2
    E = 3
    src = [1,1,2]
    tgt = [1,1,2]
  end
  julia> unique_by(g, :E, [:src, :tgt])
  ``` 
  """
  function unique_by(acset, table::Symbol, columns::Vector{Symbol})
    acset_copy = copy(acset)
    unique_by!(acset_copy, table, columns)
  end 

  """
      function compose_decapodes(decapodes_vars::Vector{Tuple{NamedDecapode{
        Any, Any, Symbol}, Vector{Symbol}}}, relation::RelationDiagram)

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
  function compose_decapodes(decapodes_vars::Vector{Tuple{NamedDecapode{
    Any, Any, Symbol}, Vector{Symbol}}}, relation::RelationDiagram)
    # TODO: Calling this function "compose" is something of a misnomer.  e.g.
    # You could very well give a vector containing a single decapode. Then the
    # effect of applying the relation would be something like renaming
    # variables. So this function should really be called oapply_decapodes (or
    # just oapply if you are willing to import oapply from
    # Catlab.WiringDiagrams).
    r = relation
    copies = decapodes_vars .|> first .|> copy
    # TODO: We should also check that the number of variables given in the
    # relation is the same as the number of symbols in the corresponding vector
    # of Vars.

    # Step -4: Determine the mapping of global ports to local ports.
    # In a RelationDiagram, this is baked into the order of rows in the Port
    # table.
    # (i.e. for each (box,junction) pair, determine whether this corresponds to
    # the 1st Var of that box, or the 2nd Var of that box, etc.)
    # tuples (box, indices of rows with that box):
    box_rows = [(b, incident(r, b, :box)) for b in unique(r[:box])]
    local_ports = [map(y -> y-idxs[begin]+1, idxs) for (b, idxs) in box_rows]
    # This is a column that one could hcat to the Ports table.
    local_ports = local_ports |> flatten |> collect

    # Step -3: Check that the types of all Vars connected by the same
    # junction are the same.
    foreach(parts(r, :Junction)) do j
      # Check that all types are equal to the first type found.
      j_idxs = incident(r, j, :junction)
      first_box = copies[r[:box][j_idxs]] |> first
      first_lp_idx = local_ports[j_idxs] |> first
      first_type = first_box[:type][first_lp_idx]
      foreach(r[:box][j_idxs], local_ports[j_idxs]) do b_idx, lp_idx
        # Get the index of the row with this name in the box's Var table.
        name = decapodes_vars[b_idx][2][lp_idx]
        local_name_idx = incident(copies[b_idx], name, [:name]) |> only
        # TODO: The `infer` type will never be here once the decapodes parsing
        # is finished being refactored. So we don't check for it here.
        # Note: This only returns the first type error found.
        copies[b_idx][:type][local_name_idx] == first_type ||
          error("The type of $(copies[b_idx][:name][local_name_idx]),
            $(copies[b_idx][:type][local_name_idx]), in decapode
            \"$(r[:name][b_idx])\" does not match the type of
            $(first_box[:name][local_name_idx]), $(first_type), in decapode
            $(r[:name][first(j_idxs)]). (Also, check that the order of the
            decapodes you supplied matches the the order you specified in the
            relation.)")
      end
    end

    # Step -2: Do namespacing.
    # Append each Var name with the name @relation gave the decapode.
    for (copy_idx, box_name) in enumerate(r[:name])
      for (name_idx, var_name) in enumerate(copies[copy_idx][:name])
        copies[copy_idx][:name][name_idx] = Symbol(
          string(box_name)*"_"*string(copies[copy_idx][:name][name_idx]))
      end
    end

    # Step -1: Write over the name fields to be what was specified by
    # @relation. (oapply cannot combine objects whose attrtypes are not equal.)
    foreach(r[:box], r[:junction], local_ports) do b_idx, j_idx, lp_idx
      # Get the index of the row with this name in the Var.
      name = decapodes_vars[b_idx][2][lp_idx]
      name = Symbol(string(r[:name][b_idx])*"_"*string(name))
      #local_name_idx = incident(copies[b_idx], name, :name) |> only
      local_name_idx = findfirst(==(name), copies[b_idx][:name]) |> only
      copies[b_idx][:name][local_name_idx] = r[:variable][j_idx]
    end

    # Step 2: Start composing.
    OpenNamedDecapodeOb, OpenNamedDecapode = OpenACSetTypes(NamedDecapode, :Var)
    
    OpenNamedDecapodes = map(copies, decapodes_vars) do curr_copy, d_vs
      variables = last(d_vs)
      FinFunctions = map(variables) do var
        FinFunction(incident(curr_copy, var, :name), nparts(curr_copy, :Var))
      end
      # TODO: It would be less brittle to pass the type information from the
      # decapodes.
      OpenNamedDecapode{Any, Any, Symbol}(curr_copy, FinFunctions...)
    end
    dec = oapply(relation, OpenNamedDecapodes) |> apex

    # Step 3: De-duplicate Op1s.
    unique_by!(dec, :Op1, [:src, :tgt, :op1])

    # Step 4: De-duplicate Op2s.
    unique_by!(dec, :Op2, [:proj1, :proj2, :res, :op2])

    # TODO: In case we have to de-duplicate tvars, this is the code.
    ## Step 5: De-duplicate TVars.
    #unique_by!(dec, :TVar, [:incl])

    # TODO: There must be a final step where you remove namespacing.

    return dec
  end

  # TODO: Is this the best way to get multiple dispatch to work?
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

  # Test that Op2s are properly de-duplicated.
  AdvectionExprBody =  quote
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
