module Examples

using CombinatorialSpaces
using Catlab
using Catlab.Present
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.CategoricalAlgebra
using Decapods.Simulations
using CombinatorialSpaces: ∧

using LinearAlgebra

import Catlab.Programs: @program
import Catlab.Graphics: to_graphviz

export dual, DWD_Funcs, sym2func, @program, to_graphviz

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
end
