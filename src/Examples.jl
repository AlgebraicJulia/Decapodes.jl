module Examples

using CombinatorialSpaces
using Catlab
using Catlab.Present
using Catlab.Programs
using Catlab.WiringDiagrams
using Catlab.Graphics
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Decapods.Simulations

import Catlab.Programs: @program
import Catlab.Graphics: to_graphviz

export dual, DWD_Funcs, sym2func, @program, to_graphviz

function dual(s::EmbeddedDeltaSet2D{O, P}) where {O, P}
  sd = EmbeddedDeltaDualComplex2D{O, eltype(P), P}(s)
  subdivide_duals!(sd, Barycenter())
  sd
end

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

sym2func(sd, k, ν) = begin
  ∂1_inds = findall(x -> x != 0, boundary(Val{2},sd) * fill(1,ntriangles(sd)))
  ∂0_inds = unique(vcat(sd[∂1_inds,:src],sd[∂1_inds,:tgt]))
  ∂1_mask = ones(Int64, ne(sd))
  ∂1_mask[∂1_inds] .= 0
  ∂0_mask = ones(Int64, nv(sd))
  ∂0_mask[∂0_inds] .= 0
  Dict(:d₀=>Dict(:operator => d(Val{0}, sd), :type => MatrixFunc()),
       :d₁=>Dict(:operator => d(Val{1}, sd), :type => MatrixFunc()),
       :dual_d₀=>Dict(:operator => dual_derivative(Val{0}, sd), :type => MatrixFunc()),
       :dual_d₁=>Dict(:operator => dual_derivative(Val{1}, sd), :type => MatrixFunc()),
       :⋆₀=>Dict(:operator => hodge_star(Val{0}, sd), :type => MatrixFunc()),
       :⋆₁=>Dict(:operator => hodge_star(Val{1}, sd), :type => MatrixFunc()),
       :⋆₂=>Dict(:operator => hodge_star(Val{2}, sd), :type => MatrixFunc()),
       :⋆₀⁻¹=>Dict(:operator => inv_hodge_star(Val{0}, sd), :type => MatrixFunc()),
       :⋆₁⁻¹=>Dict(:operator => inv_hodge_star(Val{1}, sd; hodge=DiagonalHodge()), :type => MatrixFunc()),
       :⋆₂⁻¹=>Dict(:operator => inv_hodge_star(Val{2}, sd), :type => MatrixFunc()),
       :plus => Dict(:operator => (x,y)->(x+y), :type => ElementwiseFunc()),
       :k => Dict(:operator => (x)->(k*x), :type => ElementwiseFunc()),
       :ν => Dict(:operator => (x)->(ν*x), :type => ElementwiseFunc()),
       :L₀ => Dict(:operator => (v,u)->(lie_derivative_flat(Val{2}, sd, v, u)), :type => ArbitraryFunc()),
       :Δ=>Dict(:operator => Δ(Val{0}, sd), :type => MatrixFunc()),
       :∂₁ => Dict(:operator => ∂1_mask, :type=>ConstantFunc()),
       :∂₀ => Dict(:operator => (x, y) -> x .= (∂0_mask .* y), :type=>InPlaceFunc()))
end
end
