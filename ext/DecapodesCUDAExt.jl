module DecapodesCUDAExt
using CombinatorialSpaces
using LinearAlgebra
using Base.Iterators
using Krylov
using CUDA
using CUDA.CUSPARSE
using MLStyle
import Decapodes: default_dec_cu_matrix_generate

function default_dec_cu_matrix_generate(sd, my_symbol, hodge)
  op = @match my_symbol begin

    # Regular Hodge Stars
    :⋆₀ => dec_cu_mat_hodge(0, sd, hodge)
    :⋆₁ => dec_cu_mat_hodge(Val{1}, sd, hodge)
    :⋆₂ => dec_cu_mat_hodge(2, sd, hodge)

    # Inverse Hodge Stars
    :⋆₀⁻¹ => dec_cu_mat_inverse_hodge(0, sd, hodge)
    :⋆₁⁻¹ => dec_cu_pair_inv_hodge(Val{1}, sd, hodge) # Special since Geo is a solver
    :⋆₂⁻¹ => dec_cu_mat_inverse_hodge(1, sd, hodge)

    # Differentials
    :d₀ => dec_cu_mat_differential(0, sd)
    :d₁ => dec_cu_mat_differential(1, sd)

    # Dual Differentials
    :dual_d₀ || :d̃₀ => dec_cu_mat_dual_differential(0, sd)
    :dual_d₁ || :d̃₁ => dec_cu_mat_dual_differential(1, sd)

    # Wedge Products
    :∧₀₁ => dec_cu_pair_wedge_product(Tuple{0,1}, sd)
    :∧₁₀ => dec_cu_pair_wedge_product(Tuple{1,0}, sd)
    :∧₀₂ => dec_cu_pair_wedge_product(Tuple{0,2}, sd)
    :∧₂₀ => dec_cu_pair_wedge_product(Tuple{2,0}, sd)
    :∧₁₁ => dec_cu_pair_wedge_product(Tuple{1,1}, sd)

    _ => error("Unmatched operator $my_symbol")
  end

  return op
end

# TODO: Update this to better cast hodges
function dec_cu_mat_hodge(k, sd::HasDeltaSet, hodge)
  hodge = dec_hodge_star(k, sd, hodge, Val{:CUDA})
  return (hodge, x -> hodge * x)
end

function dec_cu_mat_hodge(::Type{Val{1}}, sd::HasDeltaSet, hodge::DiagonalHodge)
  hodge = dec_hodge_star(1, sd, hodge, Val{:CUDA})
  return (hodge, x -> hodge * x)
end

function dec_cu_mat_hodge(::Type{Val{1}}, sd::HasDeltaSet, hodge::GeometricHodge)
  hodge = dec_hodge_star(Val{1}, sd, hodge, Val{:CUDA})
  return (hodge, x -> hodge * x)
end

function dec_cu_mat_inverse_hodge(k::Int, sd::HasDeltaSet, hodge)
  invhodge = dec_inv_hodge_star(k, sd, hodge, Val{:CUDA})
  return (invhodge, x -> invhodge * x)
end

"""dec_cu_pair_inv_hodge(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, _p} where _p,  ::GeometricHodge ) where float_type

This function uses Krylov.jl to compute the Geomtric Hodge Inverse for 1-Forms using GMRES. It creates an
inplace GMRES solver and returns an in-place function in the first argument and the out-of-place in the second.
"""
function dec_cu_pair_inv_hodge(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, _p} where _p,  ::GeometricHodge ) where float_type
  hdg = -1 * dec_hodge_star(Val{1}, sd, GeometricHodge(), Val{:CUDA})
  gmres_solver = GmresSolver(size(hdg, 1), size(hdg, 2), 10, CuVector{float_type})

  ((y, b) -> (y .= gmres!(gmres_solver, hdg, b).x), b -> gmres!(gmres_solver, hdg, b).x)
end

function dec_cu_pair_inv_hodge(::Type{Val{1}}, sd::HasDeltaSet, ::DiagonalHodge)
  inv_hdg = dec_inv_hodge_star(1, sd, DiagonalHodge(), Val{:CUDA})
  ((y, x) -> mul!(y, inv_hdg, x), x -> inv_hdg * x)
end

function dec_cu_mat_differential(k::Int, sd::HasDeltaSet)
  diff = dec_differential(k, sd, Val{:CUDA})
  return (diff, x -> diff * x)
end

function dec_cu_mat_dual_differential(k::Int, sd::HasDeltaSet)
  dualdiff = dec_dual_derivative(k, sd, Val{:CUDA})
  return (dualdiff, x -> dualdiff * x)
end

function dec_cu_pair_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where {k}
  val_pack = dec_p_wedge_product(Tuple{0,k}, sd, Val{:CUDA})
  ((y, α, g) -> dec_c_wedge_product!(Tuple{0,k}, y, g, α, val_pack, Val{:CUDA}),
    (α, g) -> dec_c_wedge_product(Tuple{0,k}, g, α, val_pack, Val{:CUDA}))
end

function dec_cu_pair_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where {k}
  val_pack = dec_p_wedge_product(Tuple{0,k}, sd, Val{:CUDA})
  ((y, f, β) -> dec_c_wedge_product!(Tuple{0,k}, y, f, β, val_pack, Val{:CUDA}),
    (f, β) -> dec_c_wedge_product(Tuple{0,k}, f, β, val_pack, Val{:CUDA}))
end

function dec_cu_pair_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D)
  val_pack = dec_p_wedge_product(Tuple{1,1}, sd, Val{:CUDA})
  ((y, α, β) -> dec_c_wedge_product!(Tuple{1,1}, y, α, β, val_pack, Val{:CUDA}),
    (α, β) -> dec_c_wedge_product(Tuple{1,1}, α, β, val_pack, Val{:CUDA}))
end

function dec_pair_wedge_product(::Type{Tuple{0,0}}, sd::HasDeltaSet)
  error("Replace me in compiled code with element-wise multiplication (.*)")
end

# TODO: These need to be converted into CuArrays/kernels
function dec_cu_sharp_p(sd::HasDeltaSet2D)
  ♯_m = ♯_mat(sd, AltPPSharp())
  x -> ♯_m * x
end

function dec_cu_sharp_d(sd::HasDeltaSet2D)
  ♯_m = ♯_mat(sd, LLSDDSharp())
  x -> ♯_m * x
end

function dec_cu_flat(sd::HasDeltaSet2D)
  ♭_m = ♭_mat(sd)
  x -> ♭_m * x
end
end
