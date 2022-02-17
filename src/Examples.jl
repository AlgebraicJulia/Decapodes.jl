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

import Catlab.Graphics: to_graphviz

export dual, sym2func, to_graphviz, gen_dec_rules, contract_matrices, expand_dwd, zip_dwd

function dual(s::EmbeddedDeltaSet2D{O, P}) where {O, P}
  sd = EmbeddedDeltaDualComplex2D{O, eltype(P), P}(s)
  subdivide_duals!(sd, Barycenter())
  sd
end

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
    rem_boxes_boot!(dwd, findall(b -> dwd.diagram[b, :value] == :remove, 1:nboxes(dwd)))
  end
end

""" rem_boxes_boot!

This function removes boxes from a DWD in a way that preserves port ordering.
This is only a temporary fix, and is dependent on [this
issue](https://github.com/AlgebraicJulia/Catlab.jl/issues/530) in the upstream
dependency Catlab
"""
function rem_boxes_boot!(dwd, box_to_rem)
  if isempty(box_to_rem)
    return
  end
  boxes = collect(parts(dwd.diagram, :Box))
  filter!(b -> !(b ∈ box_to_rem), boxes)
  append!(boxes, box_to_rem)
  in_ports = incident(dwd.diagram, boxes, :in_port_box)
  out_ports = incident(dwd.diagram, boxes, :out_port_box)

  in_port_trans = fill(1, nparts(dwd.diagram, :InPort))
  out_port_trans = fill(1, nparts(dwd.diagram, :OutPort))
  in_port_trans[vcat(in_ports...)] .= parts(dwd.diagram, :InPort)
  out_port_trans[vcat(out_ports...)] .= parts(dwd.diagram, :OutPort)

  set_subparts!(dwd.diagram, parts(dwd.diagram, :InPort),
                in_port_box = copy(dwd.diagram[vcat(in_ports...), :in_port_box]),
                in_port_type= copy(dwd.diagram[vcat(in_ports...), :in_port_type]))
  set_subparts!(dwd.diagram, parts(dwd.diagram, :OutPort),
                out_port_box = copy(dwd.diagram[vcat(out_ports...), :out_port_box]),
                out_port_type= copy(dwd.diagram[vcat(out_ports...), :out_port_type]))

  set_subparts!(dwd.diagram, parts(dwd.diagram, :Wire);
                src = out_port_trans[dwd.diagram[:src]],
                tgt = in_port_trans[dwd.diagram[:tgt]])
  set_subparts!(dwd.diagram, parts(dwd.diagram, :InWire),
                in_tgt = in_port_trans[dwd.diagram[:in_tgt]])
  set_subparts!(dwd.diagram, parts(dwd.diagram, :OutWire),
                out_src = out_port_trans[dwd.diagram[:out_src]])

  rem_boxes!(dwd, box_to_rem)
end

function zip_dwd!(dwd)
  updated = true
  while updated
    updated = false
    box_inputs = map(parts(dwd.diagram, :Box)) do b
      [dwd.diagram[b, :value], [(w.source.box, w.source.port) for w in in_wires(dwd, b)]]
    end
    for b1 in 1:length(box_inputs)
      for b2 in (b1+1):length(box_inputs)
        if box_inputs[b1] == box_inputs[b2]
          b1_oports = incident(dwd.diagram, b1, :out_port_box)
          b2_oports = incident(dwd.diagram, b2, :out_port_box)
          for p in 1:length(b1_oports)
            b2wires = incident(dwd.diagram, b2_oports[p], :src)
            set_subparts!(dwd.diagram, b2wires, src=b1_oports[p])
            b2outwires = incident(dwd.diagram, b2_oports[p], :out_src)
            set_subparts!(dwd.diagram, b2outwires, out_src=b1_oports[p])
          end
          rem_boxes_boot!(dwd, [b2])
          updated = true
          break
        end
      end
      if updated
        break
      end
    end
  end
end


function boundary_inds(::Type{Val{1}}, s)
  collect(findall(x -> x != 0, boundary(Val{2},s) * fill(1,ntriangles(s))))
end

function boundary_inds(::Type{Val{0}}, s)
  ∂1_inds = boundary_inds(Val{1}, s)
  unique(vcat(s[∂1_inds,:∂v0],s[∂1_inds,:∂v1]))
end

function boundary_inds(::Type{Val{2}}, s)
  ∂1_inds = boundary_inds(Val{1}, s)
  inds = map([:∂e0, :∂e1, :∂e2]) do esym
    vcat(incident(s, ∂1_inds, esym)...)
  end
  unique(vcat(inds...))
end


function bound_edges(s, ∂₀)
  te = vcat(incident(s, ∂₀, :∂v1)...)
  se = vcat(incident(s, ∂₀, :∂v0)...)
  intersect(te, se)
end

function adj_edges(s, ∂₀)
  te = vcat(incident(s, ∂₀, :∂v1)...)
  se = vcat(incident(s, ∂₀, :∂v0)...)
  unique(vcat(te, se))
end

function gen_dec_rules()
  @present ExtendedOperators(FreeExtCalc2D) begin
    X::Space
    neg₁̃::Hom(DualForm1(X), DualForm1(X)) # negative
    half₁̃::Hom(DualForm1(X), DualForm1(X)) # half
    proj₁_¹⁰₁::Hom(Form1(X) ⊗ Form0(X), Form1(X))
    proj₂_¹⁰₀::Hom(Form1(X) ⊗ Form0(X), Form0(X))
    proj₁_¹²₁::Hom(Form1(X) ⊗ Form2(X), Form1(X))
    proj₂_¹²₂::Hom(Form1(X) ⊗ Form2(X), Form2(X))
    proj₁_¹²̃₁::Hom(Form1(X) ⊗ DualForm2(X), Form1(X))
    proj₂_¹²̃₂̃::Hom(Form1(X) ⊗ DualForm2(X), DualForm2(X))
    proj₁_¹¹̃₁::Hom(Form1(X) ⊗ DualForm1(X), Form1(X))
    proj₂_¹¹̃₁̃::Hom(Form1(X) ⊗ DualForm1(X), DualForm1(X))
    proj₁_¹̃¹̃₁̃::Hom(DualForm1(X) ⊗ DualForm1(X), DualForm1(X))
    proj₂_¹̃¹̃₁̃::Hom(DualForm1(X) ⊗ DualForm1(X), DualForm1(X))
    proj₁_¹¹₁::Hom(Form1(X)⊗Form1(X), Form1(X))
    proj₂_¹¹₁::Hom(Form1(X)⊗Form1(X), Form1(X))
    sum₁̃::Hom(DualForm1(X)⊗DualForm1(X), DualForm1(X))
    sum₁::Hom(Form1(X)⊗Form1(X), Form1(X))
    L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
    L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
    i₀::Hom(Form1(X)⊗DualForm2(X), DualForm1(X))
    i₁::Hom(Form1(X)⊗DualForm1(X), DualForm0(X))
  end

  Lie0Imp = @decapode ExtendedOperators begin
    (dḞ2, dF2)::DualForm2{X}
    F1::Form1{X}

    dḞ2 ==  dual_d₁{X}(i₀(F1, dF2))
  end
  lie0_imp = diag2dwd(Lie0Imp, in_vars = [:F1, :dF2], out_vars = [:dḞ2])

  Lie1Imp = @decapode ExtendedOperators begin
    (dḞ1, dF1)::DualForm1{X}
    F1::Form1{X}
    dḞ1 == i₀(F1, dual_d₁{X}(dF1)) + dual_d₀{X}(i₁(F1, dF1))
  end
  lie1_imp = diag2dwd(Lie1Imp, in_vars = [:F1, :dF1], out_vars = [:dḞ1])

  I0Imp = @decapode ExtendedOperators begin
    F1::Form1{X}
    dF1::DualForm1{X}
    dF2::DualForm2{X}
    dF1 == ⋆₁{X}(∧₁₀{X}(F1, ⋆₀⁻¹{X}(dF2)))
  end
  i0_imp = diag2dwd(I0Imp, in_vars = [:F1, :dF2], out_vars = [:dF1])

  I1Imp = @decapode ExtendedOperators begin
    F1::Form1{X}
    dF1::DualForm1{X}
    dF0::DualForm0{X}
    dF0 == ⋆₂{X}(∧₁₁{X}(F1, ⋆₁⁻¹{X}(dF1)))
  end
  i1_imp = diag2dwd(I1Imp, in_vars = [:F1, :dF1], out_vars = [:dF0])

  δ₁Imp = @decapode ExtendedOperators begin
    F1::Form1{X}
    F0::Form0{X}
    F0 == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(F1)))
  end
  δ₁_imp = diag2dwd(δ₁Imp, in_vars = [:F1], out_vars = [:F0])

  δ₂Imp = @decapode ExtendedOperators begin
    F1::Form1{X}
    F2::Form2{X}
    F1 == ⋆₁⁻¹{X}(dual_d₀{X}(⋆₂{X}(F2)))
  end
  δ₂_imp = diag2dwd(δ₂Imp, in_vars = [:F2], out_vars = [:F1])

  Δ0Imp = @decapode ExtendedOperators begin
    (F0, Ḟ0)::Form0{X}
    Ḟ0 == δ₁{X}(d₀{X}(F0))
  end
  Δ0_imp = diag2dwd(Δ0Imp, in_vars = [:F0], out_vars = [:Ḟ0])

  Δ1Imp = @decapode ExtendedOperators begin
    (F1, Ḟ1)::Form1{X}
    Ḟ1 == δ₂{X}(d₁{X}(F1)) + d₀{X}(δ₁{X}(F1))
  end
  Δ1_imp = diag2dwd(Δ1Imp, in_vars = [:F1], out_vars = [:Ḟ1])

  Dict(:L₀ => lie0_imp, :i₀ => i0_imp, :L₁ => lie1_imp, :i₁ => i1_imp,
       :δ₁ => δ₁_imp, :δ₂ => δ₂_imp, :Δ₀ => Δ0_imp, :Δ₁ => Δ1_imp)
end
end
