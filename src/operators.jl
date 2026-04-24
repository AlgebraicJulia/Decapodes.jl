using Base.Iterators
using CombinatorialSpaces
import CombinatorialSpaces.DiscreteExteriorCalculus: DiscreteHodge
using Krylov
using LinearAlgebra
using SparseArrays

# Define mappings for default DEC operations that are not optimizable.
# --------------------------------------------------------------------
function default_dec_generate(fs::PrimalGeometricMapSeries, my_symbol::Symbol, hodge::DiscreteHodge)
  op = @match my_symbol begin
    # Inverse Laplacians
    :Δ₀⁻¹ => dec_Δ⁻¹(Val(0), fs)
    _ => default_dec_generate(finest_mesh(fs), my_symbol, hodge)
  end
end

function default_dec_generate(sd::HasDeltaSet, my_symbol::Symbol, hodge::DiscreteHodge=GeometricHodge())

  op = @match my_symbol begin

    # Misc.
    :plus => (+)
    :(-) || :neg => x -> -1 .* x
    :ln => (x -> log.(x))
    :mag => x -> norm.(x)
    :norm => x -> norm.(x)

    # Musical Isomorphisms
    :♯ᵖᵈ => dec_♯_pd(sd)
    :♯ᵖᵖ => dec_♯_pp(sd)
    :♯ᵈᵈ => dec_♯_dd(sd)
    :♭ᵈᵖ => dec_♭(sd)

    # Primal-Dual Wedge Products
    :∧ᵖᵈ₁₁ => dec_wedge_product_pd(Val(1), Val(1), sd)
    :∧ᵖᵈ₀₁ => dec_wedge_product_pd(Val(0), Val(1), sd)
    :∧ᵈᵖ₁₁ => dec_wedge_product_dp(Val(1), Val(1), sd)
    :∧ᵈᵖ₁₀ => dec_wedge_product_dp(Val(1), Val(0), sd)

    # Dual-Dual Wedge Products
    :∧ᵈᵈ₁₁ => dec_wedge_product_dd(Val(1), Val(1), sd)
    :∧ᵈᵈ₁₀ => dec_wedge_product_dd(Val(1), Val(0), sd)
    :∧ᵈᵈ₀₁ => dec_wedge_product_dd(Val(0), Val(1), sd)

    # Dual-Dual Interior Products
    :ι₁₁ => interior_product_dd(Val(1), Val(1), sd)
    :ι₁₂ => interior_product_dd(Val(1), Val(2), sd)

    # Dual-Dual Lie Derivatives
    :ℒ₁ => ℒ_dd(Val(1), Val(1), sd)

    # Dual Laplacians
    :Δᵈ₀ => Δᵈ(Val(0), sd)
    :Δᵈ₁ => Δᵈ(Val(1), sd)

    # Inverse Laplacians
    :Δ₀⁻¹ => dec_inv_lap_solver(Val(0), sd)

    _ => error("Unmatched operator $my_symbol")
  end

  return (args...) -> op(args...)
end

function default_dec_cu_generate() end;

# Define mappings for default DEC operations that are optimizable.
# ----------------------------------------------------------------
function default_dec_cu_matrix_generate() end;

function default_dec_matrix_generate(fs::PrimalGeometricMapSeries, my_symbol::Symbol, hodge::DiscreteHodge)
  op = @match my_symbol begin
    _ => default_dec_matrix_generate(finest_mesh(fs), my_symbol, hodge)
  end
end

function default_dec_matrix_generate(sd::HasDeltaSet, my_symbol::Symbol, hodge::DiscreteHodge)

  matmul(m) = (m, x -> m * x)

  op = @match my_symbol begin

    # Regular Hodge Stars
    :⋆₀ => dec_hodge_star(0, sd, hodge=hodge) |> matmul
    :⋆₁ => dec_hodge_star(1, sd, hodge=hodge) |> matmul
    :⋆₂ => dec_hodge_star(2, sd, hodge=hodge) |> matmul

    # Inverse Hodge Stars
    :⋆₀⁻¹ => dec_inv_hodge_star(0, sd, hodge) |> matmul
    :⋆₁⁻¹ => dec_pair_inv_hodge(Val(1), sd, hodge) # Special since Geo is a solver
    :⋆₂⁻¹ => dec_inv_hodge_star(1, sd, hodge) |> matmul

    # Differentials
    :d₀ => dec_differential(0, sd) |> matmul
    :d₁ => dec_differential(1, sd) |> matmul

    # Dual Differentials
    :dual_d₀ || :d̃₀ => dec_dual_derivative(0, sd) |> matmul
    :dual_d₁ || :d̃₁ => dec_dual_derivative(1, sd) |> matmul

    # Wedge Products
    :∧₀₁ => dec_pair_wedge_product(Val(0), Val(1), sd)
    :∧₁₀ => dec_pair_wedge_product(Val(1), Val(0), sd)
    :∧₀₂ => dec_pair_wedge_product(Val(0), Val(2), sd)
    :∧₂₀ => dec_pair_wedge_product(Val(2), Val(0), sd)
    :∧₁₁ => dec_pair_wedge_product(Val(1), Val(1), sd)

    :♭♯ => ♭♯_mat(sd) |> matmul

    # Averaging Operator
    :avg₀₁ => avg₀₁_mat(sd) |> matmul

     _ => error("Unmatched operator $my_symbol")
  end

  return op
end

# Special case for inverse hodge for DualForm1 to Form1
function dec_pair_inv_hodge(::Val{1}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge)
  inv_hdg = LinearAlgebra.factorize(-1 * dec_hodge_star(1, sd, GeometricHodge()))
  ((y, x) -> ldiv!(y, inv_hdg, x), x -> inv_hdg \ x)
end

function dec_pair_inv_hodge(::Val{1}, sd::HasDeltaSet, ::DiagonalHodge)
  inv_hdg = dec_inv_hodge_star(1, sd, DiagonalHodge())
  ((y, x) -> mul!(y, inv_hdg, x), x -> inv_hdg * x)
end

function dec_pair_wedge_product(::Val{k}, ::Val{0}, sd::HasDeltaSet) where {k}
  val_pack = cache_wedge(Val(0), Val(k), sd, Val(:CPU))
  ((y, α, g) -> dec_c_wedge_product!(Val(0), Val(k), y, g, α, val_pack[1], val_pack[2]),
    (α, g) -> dec_c_wedge_product(Val(0), Val(k), g, α, val_pack))
end

function dec_pair_wedge_product(::Val{0}, ::Val{k}, sd::HasDeltaSet) where {k}
  val_pack = cache_wedge(Val(0), Val(k), sd, Val(:CPU))
  ((y, f, β) -> dec_c_wedge_product!(Val(0), Val(k), y, f, β, val_pack[1], val_pack[2]),
    (f, β) -> dec_c_wedge_product(Val(0), Val(k), f, β, val_pack))
end

function dec_pair_wedge_product(::Val{1}, ::Val{1}, sd::HasDeltaSet2D)
  val_pack = cache_wedge(Val(1), Val(1), sd, Val(:CPU))
  ((y, α, β) -> dec_c_wedge_product!(Val(1), Val(1), y, α, β, val_pack[1], val_pack[2]),
    (α, β) -> dec_c_wedge_product(Val(1), Val(1), α, β, val_pack))
end

function dec_pair_wedge_product(::Val{0}, ::Val{0}, sd::HasDeltaSet)
  error("Replace me in compiled code with element-wise multiplication (.*)")
end

function dec_♯_pd(sd::HasDeltaSet1D)
  x -> ♯(sd, x, PDSharp())
end

function dec_♯_pp(sd::HasDeltaSet1D)
  x -> ♯(sd, x, PPSharp())
end

function dec_♯_pd(sd::HasDeltaSet2D)
  error("Primal-dual sharp is not yet implemented for 2D complexes.")
end

function dec_♯_pp(sd::HasDeltaSet2D)
  ♯_m = ♯_mat(sd, AltPPSharp())
  x -> ♯_m * x
end

function dec_♯_dd(sd::HasDeltaSet2D)
  ♯_m = ♯_mat(sd, LLSDDSharp())
  x -> ♯_m * x
end

function dec_♭(sd::HasDeltaSet2D)
  ♭_m = ♭_mat(sd)
  x -> ♭_m * x
end

function dec_inv_lap_solver(::Val{0}, sd::HasDeltaSet)
  inv_lap = LinearAlgebra.factorize(∇²(0, sd))
  x -> inv_lap \ x
end

function open_operators(d::SummationDecapode; kwargs...)
  e = deepcopy(d)
  open_operators!(e, kwargs...)
  return e
end

open_operators!(d::SummationDecapode; dimension::Int = 2) =
  rewrite!(d; dimension = dimension)

