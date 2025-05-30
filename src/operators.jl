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
    :Δ₀⁻¹ => dec_Δ⁻¹(Val{0}, fs)
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
    :Δᵈ₀ => Δᵈ(Val{0}, sd)
    :Δᵈ₁ => Δᵈ(Val{1}, sd)

    # Inverse Laplacians
    :Δ₀⁻¹ => dec_inv_lap_solver(Val{0}, sd)

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
    :⋆₁⁻¹ => dec_pair_inv_hodge(Val{1}, sd, hodge) # Special since Geo is a solver
    :⋆₂⁻¹ => dec_inv_hodge_star(1, sd, hodge) |> matmul

    # Differentials
    :d₀ => dec_differential(0, sd) |> matmul
    :d₁ => dec_differential(1, sd) |> matmul

    # Dual Differentials
    :dual_d₀ || :d̃₀ => dec_dual_derivative(0, sd) |> matmul
    :dual_d₁ || :d̃₁ => dec_dual_derivative(1, sd) |> matmul

    # Wedge Products
    :∧₀₁ => dec_pair_wedge_product(Tuple{0,1}, sd)
    :∧₁₀ => dec_pair_wedge_product(Tuple{1,0}, sd)
    :∧₀₂ => dec_pair_wedge_product(Tuple{0,2}, sd)
    :∧₂₀ => dec_pair_wedge_product(Tuple{2,0}, sd)
    :∧₁₁ => dec_pair_wedge_product(Tuple{1,1}, sd)

    :♭♯ => ♭♯_mat(sd) |> matmul

    # Averaging Operator
    :avg₀₁ => avg₀₁_mat(sd) |> matmul

     _ => error("Unmatched operator $my_symbol")
  end

  return op
end

# Special case for inverse hodge for DualForm1 to Form1
function dec_pair_inv_hodge(::Type{Val{1}}, sd::AbstractDeltaDualComplex2D, ::GeometricHodge)
  inv_hdg = LinearAlgebra.factorize(-1 * dec_hodge_star(1, sd, GeometricHodge()))
  ((y, x) -> ldiv!(y, inv_hdg, x), x -> inv_hdg \ x)
end

function dec_pair_inv_hodge(::Type{Val{1}}, sd::HasDeltaSet, ::DiagonalHodge)
  inv_hdg = dec_inv_hodge_star(1, sd, DiagonalHodge())
  ((y, x) -> mul!(y, inv_hdg, x), x -> inv_hdg * x)
end

function dec_pair_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet) where {k}
  val_pack = cache_wedge(Tuple{0,k}, sd, Val{:CPU})
  ((y, α, g) -> dec_c_wedge_product!(Tuple{0,k}, y, g, α, val_pack[1], val_pack[2]),
    (α, g) -> dec_c_wedge_product(Tuple{0,k}, g, α, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet) where {k}
  val_pack = cache_wedge(Tuple{0,k}, sd, Val{:CPU})
  ((y, f, β) -> dec_c_wedge_product!(Tuple{0,k}, y, f, β, val_pack[1], val_pack[2]),
    (f, β) -> dec_c_wedge_product(Tuple{0,k}, f, β, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D)
  val_pack = cache_wedge(Tuple{1,1}, sd, Val{:CPU})
  ((y, α, β) -> dec_c_wedge_product!(Tuple{1,1}, y, α, β, val_pack[1], val_pack[2]),
    (α, β) -> dec_c_wedge_product(Tuple{1,1}, α, β, val_pack))
end

function dec_pair_wedge_product(::Type{Tuple{0,0}}, sd::HasDeltaSet)
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

function dec_inv_lap_solver(::Type{Val{0}}, sd::HasDeltaSet)
  inv_lap = LinearAlgebra.factorize(∇²(0, sd))
  x -> inv_lap \ x
end

function open_operators(d::SummationDecapode; kwargs...)
  e = deepcopy(d)
  open_operators!(e, kwargs...)
  return e
end

function open_operators!(d::SummationDecapode; dimension::Int=2)
  op1_remove_stack = Vector{Int}()
  for op1_idx in parts(d, :Op1)
    op1_src = d[op1_idx, :src]
    op1_tgt = d[op1_idx, :tgt]
    op1_name = d[op1_idx, :op1]


    remove_op1 = @match (op1_name, dimension) begin
      (:Δ₀, 1) => begin
        add_De_Rham_1D!(Val{0}, d, op1_src, op1_tgt)
        true
      end
      (:Δ₁, 1) => begin
        add_De_Rham_1D!(Val{1}, d, op1_src, op1_tgt)
        true
      end
      (:Δ₀, 2) => begin
        add_De_Rham_2D!(Val{0}, d, op1_src, op1_tgt)
        true
      end
      (:Δ₁, 2) => begin
        add_De_Rham_2D!(Val{1}, d, op1_src, op1_tgt)
        true
      end
      (:Δ₂, 2) => begin
        add_De_Rham_2D!(Val{2}, d, op1_src, op1_tgt)
        true
      end
      (:δ₁, _) => begin
        add_Codiff!(d, op1_src, op1_tgt)
        true
      end
      (:δ₂, _) => begin
        add_Codiff!(d, op1_src, op1_tgt)
        true
      end
      _ => false
    end
    remove_op1 && push!(op1_remove_stack, op1_idx)
  end

  op2_remove_stack = Vector{Int}()
  for op2_idx in parts(d, :Op2)
    op2_proj1 = d[op2_idx, :proj1]
    op2_proj2 = d[op2_idx, :proj2]
    op2_res = d[op2_idx, :res]
    op2_name = d[op2_idx, :op2]

    remove_op2 = 0

    ## Make substitution for complex operator into components
    remove_op2 = @match (op2_name, dimension) begin
      (:i₁, 1) => begin
        add_Inter_Prod_1D!(Val{1}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:i₁, 2) => begin
        add_Inter_Prod_2D!(Val{1}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:i₂, 2) => begin
        add_Inter_Prod_2D!(Val{2}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:L₀, 1) => begin
        add_Lie_1D!(Val{0}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:L₁, 1) => begin
        add_Lie_1D!(Val{1}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:L₀, 2) => begin
        add_Lie_2D!(Val{0}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:L₁, 2) => begin
        add_Lie_2D!(Val{1}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      (:L₂, 2) => begin
        add_Lie_2D!(Val{2}, d, op2_proj1, op2_proj2, op2_res)
        true
      end
      _ => false
    end

    ## If sub was made, add the original operator to be removed
    remove_op2 && push!(op2_remove_stack, op2_idx)
  end

  ## Remove all subbed operators
  rem_parts!(d, :Op1, op1_remove_stack)
  rem_parts!(d, :Op2, op2_remove_stack)

  ## Add unique names for all newly added variables
  fill_names!(d, lead_symbol=Symbol("Gensim_Var_"))
end

function add_Inter_Prod!(d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
  ## Adds the hodge Dual to Primal
  inv_hodge_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=proj1_Inter, tgt=inv_hodge_tgt, op1=:⋆)

  ## Adds the wedge between Primal and Primal
  wedge_res = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op2, proj1=inv_hodge_tgt, proj2=proj2_Inter, res=wedge_res, op2=:∧)

  ## Adds the hodge Primal to Dual
  add_part!(d, :Op1, src=wedge_res, tgt=res_Inter, op1=:⋆)
end

function add_Inter_Prod_1D!(::Type{Val{1}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
  add_Inter_Prod!(d, proj1_Inter, proj2_Inter, res_Inter)
end

function add_Inter_Prod_2D!(::Type{Val{1}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
  ## Takes generic interior product
  pos_inter_prod = add_part!(d, :Var, type=:infer, name=nothing)
  add_Inter_Prod!(d, proj1_Inter, proj2_Inter, pos_inter_prod)

  ## Outputs negated value
  add_part!(d, :Op1, src=pos_inter_prod, tgt=res_Inter, op1=:neg)
end

function add_Inter_Prod_2D!(::Type{Val{2}}, d::SummationDecapode, proj1_Inter::Int, proj2_Inter::Int, res_Inter::Int)
  add_Inter_Prod!(d, proj1_Inter, proj2_Inter, res_Inter)
end

function add_Lie_1D!(::Type{Val{0}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

  ## Outputs result of dual derivative Dual0 to Dual1
  dual_d_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=proj2_Lie, tgt=dual_d_tgt, op1=:d)

  ## Takes interior product of Primal1 and Dual1 to Dual0
  add_Inter_Prod_1D!(Val{1}, d, dual_d_tgt, proj1_Lie, res_Lie)
end

function add_Lie_1D!(::Type{Val{1}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

  ## Takes interior product of Primal1 and Dual1 to Dual0
  inter_product_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_Inter_Prod_1D!(Val{1}, d, proj2_Lie, proj1_Lie, inter_product_tgt)

  ## Outputs result of dual derivative Dual0 to Dual1
  add_part!(d, :Op1, src=inter_product_tgt, tgt=res_Lie, op1=:d)
end

function add_Lie_2D!(::Type{Val{0}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

  ## Outputs result of dual derivative Dual0 to Dual1
  dual_d_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=proj2_Lie, tgt=dual_d_tgt, op1=:d)

  ## Takes interior product of Primal1 and Dual1 to Dual0
  add_Inter_Prod_1D!(Val{1}, d, dual_d_tgt, proj1_Lie, res_Lie)
end

function add_Lie_2D!(::Type{Val{1}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

  ## Outputs result of dual derivative Dual1 to Dual2
  dual_d_1_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=proj2_Lie, tgt=dual_d_1_tgt, op1=:d)

  ## Takes interior product of Primal1 and Dual2 to Dual1
  inter_product_2_res = add_part!(d, :Var, type=:infer, name=nothing)
  add_Inter_Prod_2D!(Val{2}, d, dual_d_1_tgt, proj1_Lie, inter_product_2_res)

  ## Takes interior product of Primal1 and Dual1 to Dual0
  inter_product_1_res = add_part!(d, :Var, type=:infer, name=nothing)
  add_Inter_Prod_2D!(Val{1}, d, proj2_Lie, proj1_Lie, inter_product_1_res)

  ## Outputs result of dual derivative Dual0 to Dual1
  dual_d_0_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=inter_product_1_res, tgt=dual_d_0_tgt, op1=:d)

  ## Outputs sum of both dual_d_0 and inter_product_2
  summation_tgt = add_part!(d, :Σ, sum=res_Lie)

  add_part!(d, :Summand, summand=inter_product_2_res, summation=summation_tgt)
  add_part!(d, :Summand, summand=dual_d_0_tgt, summation=summation_tgt)
end

function add_Lie_2D!(::Type{Val{2}}, d::SummationDecapode, proj1_Lie::Int, proj2_Lie::Int, res_Lie::Int)

  ## Takes interior product of Primal1 and Dual2 to Dual1
  inter_product_2_tgt = add_part!(d, :Var, type=:infer, name=nothing)
  add_Inter_Prod_2D!(Val{2}, d, proj2_Lie, proj1_Lie, inter_product_2_tgt)

  ## Outputs result of dual derivative Dual1 to Dual2
  add_part!(d, :Op1, src=inter_product_2_tgt, tgt=res_Lie, op1=:d)
end

function add_Codiff!(d::SummationDecapode, src_Codiff::Int, tgt_Codiff::Int)
  hodge_star_first = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=src_Codiff, tgt=hodge_star_first, op1=:⋆)

  exterior_deriv = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=hodge_star_first, tgt=exterior_deriv, op1=:d)

  add_part!(d, :Op1, src=exterior_deriv, tgt=tgt_Codiff, op1=:⋆)
end

function add_De_Rham_1D!(::Type{Val{0}}, d::SummationDecapode, src_De_Rham::Int, tgt_De_Rham::Int)
  exterior_deriv = add_part!(d, :Var, type=:infer, name=nothing)
  add_part!(d, :Op1, src=src_De_Rham, tgt=exterior_deriv, op1=:d)

  add_Codiff!(d, exterior_deriv, tgt_De_Rham)
end

function add_De_Rham_1D!(::Type{Val{1}}, d::SummationDecapode, src_De_Rham::Int, tgt_De_Rham::Int)
  codiff = add_part!(d, :Var, type=:infer, name=nothing)
  add_Codiff!(d, src_De_Rham, codiff)

  add_part!(d, :Op1, src=codiff, tgt=tgt_De_Rham, op1=:d)
end

function add_De_Rham_2D!(::Type{Val{0}}, d::SummationDecapode, src_De_Rham::Int, tgt_De_Rham::Int)
  add_De_Rham_1D!(Val{0}, d, src_De_Rham, tgt_De_Rham)
end

function add_De_Rham_2D!(::Type{Val{1}}, d::SummationDecapode, src_De_Rham::Int, tgt_De_Rham::Int)
  sum_part_1 = add_part!(d, :Var, type=:infer, name=nothing)
  add_De_Rham_2D!(Val{0}, d, src_De_Rham, sum_part_1)

  sum_part_2 = add_part!(d, :Var, type=:infer, name=nothing)
  add_De_Rham_2D!(Val{2}, d, src_De_Rham, sum_part_2)

  summation_tgt = add_part!(d, :Σ, sum=tgt_De_Rham)

  add_part!(d, :Summand, summand=sum_part_1, summation=summation_tgt)
  add_part!(d, :Summand, summand=sum_part_2, summation=summation_tgt)
end

function add_De_Rham_2D!(::Type{Val{2}}, d::SummationDecapode, src_De_Rham::Int, tgt_De_Rham::Int)
  add_De_Rham_1D!(Val{1}, d, src_De_Rham, tgt_De_Rham)
end
