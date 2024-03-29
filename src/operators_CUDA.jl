using CombinatorialSpaces
using LinearAlgebra
using Base.Iterators
using Catlab
using CUDA
using CUDA.CUSPARSE

function default_cu_dec_matrix_generate(sd, my_symbol, hodge)
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

    # Primal-Dual Wedge Products
    :∧ᵖᵈ₁₁ => dec_wedge_product_pd(Tuple{1,1}, sd)
    :∧ᵖᵈ₀₁ => dec_wedge_product_pd(Tuple{0,1}, sd)
    :∧ᵈᵖ₁₁ => dec_wedge_product_dp(Tuple{1,1}, sd)
    :∧ᵈᵖ₁₀ => dec_wedge_product_dp(Tuple{1,0}, sd)

    # Dual-Dual Wedge Products
    :∧ᵈᵈ₁₁ => dec_wedge_product_dd(Tuple{1,1}, sd)
    :∧ᵈᵈ₁₀ => dec_wedge_product_dd(Tuple{1,0}, sd)
    :∧ᵈᵈ₀₁ => dec_wedge_product_dd(Tuple{0,1}, sd)

    # Dual-Dual Interior Products
    :ι₁₁ => interior_product_dd(Tuple{1,1}, sd)
    :ι₁₂ => interior_product_dd(Tuple{1,2}, sd)

    # Dual-Dual Lie Derivatives
    :ℒ₁ => ℒ_dd(Tuple{1,1}, sd)

    # Dual Laplacians
    :Δᵈ₀ => Δᵈ(Val{0},sd)
    :Δᵈ₁ => Δᵈ(Val{1},sd)

    # Musical Isomorphisms
    :♯ => dec_♯_p(sd)
    :♯ᵈ => dec_♯_d(sd)

    :♭ => dec_♭(sd)

    :neg => x -> -1 .* x
     _ => error("Unmatched operator $my_symbol")
  end

  return op
end

# TODO: Update this to better cast hodges
function dec_mat_hodge(k, sd::HasDeltaSet, hodge)
  hodge = dec_hodge_star(k, sd, hodge=hodge)
  if(hodge isa Diagonal)
    hodge = CuArray(hodge)
  else
    hodge = CuSparseMatrixCSC(hodge)
  end
  return (hodge, x -> hodge * x)
end

function dec_mat_inverse_hodge(k::Int, sd::HasDeltaSet, hodge)
  invhodge = CuArray(dec_inv_hodge_star(k, sd, hodge=hodge))
  return (invhodge, x -> invhodge * x)
end

# Special case for inverse hodge for DualForm1 to Form1
function dec_pair_inv_hodge(::Type{Val{1}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge)
  inv_hdg = LinearAlgebra.factorize(-1 * CuSparseMatrixCSC(dec_hodge_star(1, sd, GeometricHodge())))
  ((y, x) -> ldiv!(y, inv_hdg, x), x -> inv_hdg \ x)
end

function dec_pair_inv_hodge(::Type{Val{1}}, sd::HasDeltaSet, ::DiagonalHodge)
  inv_hdg = CuArray(dec_inv_hodge_star(1, sd, DiagonalHodge()))
  ((y, x) -> mul!(y, inv_hdg, x), x -> inv_hdg * x)
end

function dec_mat_differential(k::Int, sd::HasDeltaSet)
  diff = CuSparseMatrixCSC(dec_differential(k, sd))
  return (diff, x -> diff * x)
end

function dec_mat_dual_differential(k::Int, sd::HasDeltaSet)
  dualdiff = CuSparseMatrixCSC(dec_dual_derivative(k, sd))
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