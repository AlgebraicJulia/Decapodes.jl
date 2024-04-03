using CombinatorialSpaces
using LinearAlgebra
using Base.Iterators
using Catlab
using Krylov

function default_dec_matrix_generate(sd, my_symbol, hodge, code_target::gen)
  op = @match my_symbol begin

    # Regular Hodge Stars
    :⋆₀ => dec_mat_hodge(0, sd, hodge)
    :⋆₁ => dec_mat_hodge(1, sd, hodge)
    :⋆₂ => dec_mat_hodge(2, sd, hodge)

    # Inverse Hodge Stars
    :⋆₀⁻¹ => dec_mat_inverse_hodge(0, sd, hodge)
    :⋆₁⁻¹ => dec_pair_inv_hodge(Val{1}, sd, hodge) # Special since Geo is a solver
    :⋆₂⁻¹ => dec_mat_inverse_hodge(1, sd, hodge)

    # Differentials
    :d₀ => dec_mat_differential(0, sd)
    :d₁ => dec_mat_differential(1, sd)

    # Dual Differentials
    :dual_d₀ || :d̃₀ => dec_mat_dual_differential(0, sd)
    :dual_d₁ || :d̃₁ => dec_mat_dual_differential(1, sd)

    # Wedge Products
    :∧₀₁ => dec_pair_wedge_product(Tuple{0,1}, sd)
    :∧₁₀ => dec_pair_wedge_product(Tuple{1,0}, sd)
    :∧₀₂ => dec_pair_wedge_product(Tuple{0,2}, sd)
    :∧₂₀ => dec_pair_wedge_product(Tuple{2,0}, sd)
    :∧₁₁ => dec_pair_wedge_product(Tuple{1,1}, sd)

    _ => error("Unmatched operator $my_symbol")
  end

  return op
end

# TODO: Update this to better cast hodges
function dec_mat_hodge(k, sd::HasDeltaSet, hodge)
  hodge = dec_hodge_star(k, sd, hodge, Val{CUDA})
  return (hodge, x -> hodge * x)
end

function dec_mat_inverse_hodge(k::Int, sd::HasDeltaSet, hodge)
  invhodge = dec_inv_hodge_star(k, sd, hodge, Val{CUDA})
  return (invhodge, x -> invhodge * x)
end

# Special case for inverse hodge for DualForm1 to Form1
# TODO: This should be changed to use Krylov, need to figure out inplace version of this
function dec_pair_inv_hodge(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, _p},  ::GeometricHodge where _p) where float_type
  hdg = -1 * dec_hodge_star(1, sd, GeometricHodge(), Val{CUDA})
  # TODO: Figure out what a good number for this memory value is
  gmres_solver = GmresSolver(size(hdg, 1), size(hdg, 2), 10, CuVector{float_type})

  ((y, b) -> (y .= gmres!(gmres_solver, hdg, b).x), b -> gmres!(gmres_solver, hdg, b).x)
end

function dec_pair_inv_hodge(::Type{Val{1}}, sd::HasDeltaSet, ::DiagonalHodge)
  inv_hdg = dec_inv_hodge_star(1, sd, DiagonalHodge(), Val{CUDA})
  ((y, x) -> mul!(y, inv_hdg, x), x -> inv_hdg * x)
end

function dec_mat_differential(k::Int, sd::HasDeltaSet)
  diff = dec_differential(k, sd, Val{CUDA})
  return (diff, x -> diff * x)
end

function dec_mat_dual_differential(k::Int, sd::HasDeltaSet)
  dualdiff = dec_dual_derivative(k, sd, Val{CUDA})
  return (dualdiff, x -> dualdiff * x)
end

function dec_pair_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where {k}
  val_pack = dec_cu_p_wedge_product(Tuple{0,k}, sd)
  ((y, α, g) -> dec_cu_c_wedge_product!(Tuple{0,k}, y, g, α, val_pack),
    (α, g) -> dec_cu_c_wedge_product(Tuple{0,k}, g, α, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where {k}
  val_pack = dec_cu_p_wedge_product(Tuple{0,k}, sd)
  ((y, f, β) -> dec_cu_c_wedge_product!(Tuple{0,k}, y, f, β, val_pack),
    (f, β) -> dec_cu_c_wedge_product(Tuple{0,k}, f, β, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D)
  val_pack = dec_cu_p_wedge_product(Tuple{1,1}, sd)
  ((y, α, β) -> dec_cu_c_wedge_product!(Tuple{1,1}, y, α, β, val_pack),
    (α, β) -> dec_cu_c_wedge_product(Tuple{1,1}, α, β, val_pack))
end

# TODO: These need to be converted into CuArrays/kernels
function dec_♯_p(sd::HasDeltaSet2D)
  ♯_m = ♯_mat(sd, AltPPSharp())
  x -> ♯_m * x
end

function dec_♯_d(sd::HasDeltaSet2D)
  ♯_m = ♯_mat(sd, LLSDDSharp())
  x -> ♯_m * x
end

function dec_♭(sd::HasDeltaSet2D)
  ♭_m = ♭_mat(sd)
  x -> ♭_m * x
end