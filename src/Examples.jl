module Examples

using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Catlab
using Catlab.Present
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.CategoricalAlgebra
using Decapodes.Simulations
using Decapodes.Diagrams
using Decapodes.Schedules
using CombinatorialSpaces: ∧

using LinearAlgebra

import Catlab.Programs: @program
import Catlab.Graphics: to_graphviz

export dual, DWD_Funcs, sym2func, @program, to_graphviz, gen_dec_rules

function dual(s::EmbeddedDeltaSet2D{O, P}) where {O, P}
  sd = EmbeddedDeltaDualComplex2D{O, eltype(P), P}(s)
  subdivide_duals!(sd, Barycenter())
  sd
end

#=
@present DWD_Funcs(FreeSymmetricMonoidalCategory) begin
    Scalar::Ob
    Form0::Ob
    Form1::Ob
    Form2::Ob
    DForm0::Ob
    DForm1::Ob
    DForm2::Ob
    d0::Hom(Form0, Form1)
    d1::Hom(Form1, Form2)
    d̃0::Hom(DForm0, DForm1)
    d̃1::Hom(DForm1, DForm2)
    star0::Hom(Form0, DForm2)
    star1::Hom(Form1, DForm1)
    star_inv0::Hom(DForm2, Form0)
    star_inv1::Hom(DForm1, Form1)
    star_inv2::Hom(DForm0, Form2)
    plus0::Hom(Form0⊗Form0, Form0)
    plus1::Hom(Form1⊗Form1, Form1)
    dplus1::Hom(DForm1⊗DForm1, DForm1)
    dplus2::Hom(DForm2⊗DForm2, DForm2)
    prod0::Hom(Form0⊗Form0, Form0)
    prod1::Hom(Form1⊗Form1, Form1)
    scale0::Hom(Scalar⊗Form0, Form0)
    scale1::Hom(Scalar⊗Form1, Form1)
    ∂0::Hom(munit(), Form0)
    ∂1::Hom(munit(), Form1)
    z0::Hom(munit(), Form0)
    z1::Hom(munit(), Form1)
    z1::Hom(munit(), Form2)
    k::Hom(munit(), Scalar)
    ℒ::Hom(Form1⊗DForm2, DForm2)
    Δ::Hom(Form0, Form0)
end
=#
sym2func(sd) = begin
  Dict(:d₀=>Dict(:operator => d(Val{0}, sd), :type => MatrixFunc()),
       :d₁=>Dict(:operator => d(Val{1}, sd), :type => MatrixFunc()),
       :dual_d₀=>Dict(:operator => dual_derivative(Val{0}, sd),
                      :type => MatrixFunc()),
       :dual_d₁=>Dict(:operator => dual_derivative(Val{1}, sd),
                      :type => MatrixFunc()),
       :⋆₀=>Dict(:operator => hodge_star(Val{0}, sd), :type => MatrixFunc()),
       :⋆₁=>Dict(:operator => hodge_star(Val{1}, sd), :type => MatrixFunc()),
       :⋆₂=>Dict(:operator => hodge_star(Val{2}, sd), :type => MatrixFunc()),
       :⋆₀⁻¹=>Dict(:operator => inv_hodge_star(Val{0}, sd), :type => MatrixFunc()),
       :⋆₁⁻¹=>Dict(:operator => inv_hodge_star(Val{1}, sd; hodge=DiagonalHodge()),
                   :type => MatrixFunc()),
       :⋆₂⁻¹=>Dict(:operator => inv_hodge_star(Val{2}, sd), :type => MatrixFunc()),
       :∧₁₀=>Dict(:operator => (γ, α,β)->(γ .= ∧(Tuple{1,0}, sd, α, β)),
                  :type=>InPlaceFunc()),
       :∧₁₁=>Dict(:operator => (γ, α,β)->(γ .= ∧(Tuple{1,1}, sd, α, β)),
                  :type=>InPlaceFunc()),
       :plus => Dict(:operator => (x,y)->(x+y), :type => ElementwiseFunc()),
       :L₀ => Dict(:operator => (v,u)->(lie_derivative_flat(Val{2}, sd, v, u)),
                   :type => ArbitraryFunc()))
#       :Δ₀=>Dict(:operator => Δ(Val{0}, sd), :type => MatrixFunc()),
#       :Δ₁=>Dict(:operator => Δ(Val{1}, sd), :type => MatrixFunc()),
end

function expand_dwd(dwd_orig, patterns)
  dwd = deepcopy(dwd_orig)
  has_updated = true
  while(has_updated)
    updates = map(1:nboxes(dwd)) do b
      v = box(dwd, b).value
      if v ∈ keys(patterns)
        dwd.diagram[b, :value] = patterns[v]
        true
      else
        false
      end
    end
    dwd = substitute(dwd)
    has_updated = any(updates)
  end
  dwd
end

function contract_matrices!(dwd, s2f)
  has_updated = true
  is_matrix(b) = begin
    b >= 1 || return false
    v = Symbol(first(split("$(box(dwd, b).value)", "⋅")))
    v != :remove && (s2f[v][:type] isa MatrixFunc) &&
      length(input_ports(dwd, b)) == 1 &&
      length(output_ports(dwd, b)) == 1
  end
  while(has_updated)
    has_updated = false
    # Join parallel matrix boxes and mark for removal
    for b in 1:nboxes(dwd)
      if is_matrix(b)
        ows = out_wires(dwd, b)
        if length(ows) == 1
          ow = ows[1]
          if is_matrix(ow.target.box)
            has_updated = true
            tb = ow.target.box
            dwd.diagram[b, :value] = Symbol(dwd.diagram[tb, :value], :⋅, dwd.diagram[b, :value])
            dwd.diagram[tb, :value] = :remove
            add_wires!(dwd, map(w -> Wire(w.value, Port(b, OutputPort, 1), w.target),
                                out_wires(dwd, tb)))
            dwd.diagram[incident(dwd.diagram, b, :out_port_box)[1], :out_port_type] = output_ports(dwd, tb)[1]
          end
        end
      end
    end
    rem_boxes!(dwd, findall(b -> dwd.diagram[b, :value] == :remove, 1:nboxes(dwd)))
  end
end


function boundary_inds(::Type{Val{1}}, s)
  collect(findall(x -> x != 0, boundary(Val{2},s) * fill(1,ntriangles(s))))
end

function boundary_inds(::Type{Val{0}}, s)
  ∂1_inds = boundary_inds(Val{1}, s)
  unique(vcat(s[∂1_inds,:src],s[∂1_inds,:tgt]))
end

function boundary_inds(::Type{Val{2}}, s)
  ∂1_inds = boundary_inds(Val{1}, s)
  inds = map([:∂e0, :∂e1, :∂e2]) do esym
    vcat(incident(s, ∂1_inds, esym)...)
  end
  unique(vcat(inds...))
end


function bound_edges(s, ∂₀)
  te = vcat(incident(s, ∂₀, :tgt)...)
  se = vcat(incident(s, ∂₀, :src)...)
  intersect(te, se)
end

function adj_edges(s, ∂₀)
  te = vcat(incident(s, ∂₀, :tgt)...)
  se = vcat(incident(s, ∂₀, :src)...)
  unique(vcat(te, se))
end

function gen_dec_rules()
  @present ExtendedOperators(FreeExtCalc2D) begin
    X::Space
    F0::Hom(munit(), Form0(X))
    F1::Hom(munit(), Form1(X))
    F2::Hom(munit(), Form2(X))
    dF0::Hom(munit(), DualForm0(X))
    dF1::Hom(munit(), DualForm1(X))
    dF2::Hom(munit(), DualForm2(X))
    neg::Hom(DualForm1(X), DualForm1(X)) # negative
    half::Hom(DualForm1(X), DualForm1(X)) # half
    L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
    L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
    i₀::Hom(Form1(X)⊗DualForm2(X), DualForm1(X))
    i₁::Hom(Form1(X)⊗DualForm1(X), DualForm0(X))
  end

  @present Lie0Imp <: ExtendedOperators begin
    dF2 ⋅ ∂ₜ(DualForm2(X)) == (F1 ⊗ dF2) ⋅i₀ ⋅ dual_d₁(X)
  end

  @present Lie1Imp <: ExtendedOperators begin
    dF1 ⋅ ∂ₜ(DualForm1(X)) == (F1 ⊗ (dF1 ⋅ dual_d₁(X))) ⋅ i₀ + (F1 ⊗ dF1) ⋅ i₁ ⋅ dual_d₀(X)
  end

  @present I0Imp <: ExtendedOperators begin
    dF1 ⋅ ∂ₜ(DualForm1(X)) == (F1 ⊗ (dF2 ⋅ ⋆₀⁻¹(X))) ⋅ ∧₁₀(X) ⋅ ⋆₁(X)
  end

  @present I1Imp <: ExtendedOperators begin
    dF0 ⋅ ∂ₜ(DualForm0(X)) == (F1 ⊗ (dF1 ⋅ ⋆₁⁻¹(X))) ⋅ ∧₁₁(X) ⋅ ⋆₂(X)
  end

  @present δ₁Imp <: ExtendedOperators begin
    F0 ⋅ ∂ₜ(Form0(X)) == F1 ⋅ ⋆₁(X) ⋅ dual_d₁(X) ⋅ ⋆₀⁻¹(X)
  end

  @present δ₂Imp <: ExtendedOperators begin
    F1 ⋅ ∂ₜ(Form1(X)) == F2 ⋅ ⋆₂(X) ⋅ dual_d₀(X) ⋅ ⋆₁⁻¹(X)
  end

  @present Δ0Imp <: ExtendedOperators begin
    F0 ⋅ ∂ₜ(Form0(X)) == F0 ⋅ d₀(X) ⋅ δ₁(X)
  end

  @present Δ1Imp <: ExtendedOperators begin
    F1 ⋅ ∂ₜ(Form1(X)) == F1 ⋅ (d₁(X) ⋅ δ₂(X) + δ₁(X) ⋅ d₀(X))
  end

  # L₀

  lie0_imp_diag = eq_to_diagrams(Lie0Imp)
  lie0_imp = diag2dwd(lie0_imp_diag)
  tmp = lie0_imp.diagram[1, :outer_in_port_type]
  lie0_imp.diagram[1, :outer_in_port_type] = lie0_imp.diagram[2, :outer_in_port_type]
  lie0_imp.diagram[2, :outer_in_port_type] = tmp
  tmp = lie0_imp.diagram[1, :in_src]
  lie0_imp.diagram[1, :in_src] = lie0_imp.diagram[2, :in_src]
  lie0_imp.diagram[2, :in_src] = tmp

  # L₁

  lie1_imp_diag = eq_to_diagrams(Lie1Imp)
  lie1_imp = diag2dwd(lie1_imp_diag)
  tmp = lie1_imp.diagram[1, :outer_in_port_type]
  lie1_imp.diagram[1, :outer_in_port_type] = lie1_imp.diagram[2, :outer_in_port_type]
  lie1_imp.diagram[2, :outer_in_port_type] = tmp
  lie1_imp.diagram[1, :in_src] = 2
  lie1_imp.diagram[2, :in_src] = 1
  lie1_imp.diagram[3, :in_src] = 1
  lie1_imp.diagram[4, :in_src] = 2

  # i₀

  i0_imp_diag = eq_to_diagrams(I0Imp)
  i0_imp = diag2dwd(i0_imp_diag)
  rem_part!(i0_imp.diagram, :OuterInPort, 1)

  # i₁

  i1_imp_diag = eq_to_diagrams(I1Imp)
  i1_imp = diag2dwd(i1_imp_diag)
  rem_part!(i1_imp.diagram, :OuterInPort, 1)

  # δ₁

  δ₁_imp_diag = eq_to_diagrams(δ₁Imp)
  δ₁_imp = diag2dwd(δ₁_imp_diag)
  rem_part!(δ₁_imp.diagram, :OuterInPort, 1)

  # δ₂

  δ₂_imp_diag = eq_to_diagrams(δ₂Imp)
  δ₂_imp = diag2dwd(δ₂_imp_diag)
  rem_part!(δ₂_imp.diagram, :OuterInPort, 1)

  # Δ₀

  Δ0_imp_diag = eq_to_diagrams(Δ0Imp)
  Δ0_imp = diag2dwd(Δ0_imp_diag)

  # Δ₁

  Δ1_imp_diag = eq_to_diagrams(Δ1Imp)
  Δ1_imp = diag2dwd(Δ1_imp_diag)

  Dict(:L₀ => lie0_imp, :i₀ => i0_imp, :L₁ => lie1_imp, :i₁ => i1_imp,
       :δ₁ => δ₁_imp, :δ₂ => δ₂_imp, :Δ₀ => Δ0_imp, :Δ₁ => Δ1_imp)
end
end
