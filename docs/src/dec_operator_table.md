# Discrete Exterior Calculus Operator Table

Preface: This document is meant to provide simplified descriptions of the operators of Discrete Exterior Calculus. The reader can find more rigorous definitions in the referenced sources for each operator.

## Boundary

### Signature

Decapodes.jl

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ∂(s::HasDeltaSet, x::DualChain{n}) where n
    dual_boundary(n::Int, s::HasDeltaSet, args...)
    dual_boundary(::Type{Val{n}}, s::HasDeltaSet, args...) where n
    dual_boundary_nz(::Type{Val{1}}, s::AbstractDeltaDualComplex1D, x::Int)
    dual_boundary_nz(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, x::Int)
    dual_boundary_nz(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, x::Int)

CombinatorialSpaces.jl, FastDEC.jl

    dec_boundary(n::Int, sd::HasDeltaSet)
    dec_p_boundary(::Type{Val{k}}, sd::HasDeltaSet; negate::Bool=false) where {k}

### Description

The boundary operator acts on the geometry of the mesh, linearly mapping k-chains to (k-1)-chains. It relates to the important exterior derivative operator via duality, and helps compute the operator through Generalized Stokes' Theorem.

This operator inputs k-chains and outputs (k-1)-chains (or k-manifolds to (k-1)-manifolds in the smooth case). For example, it would input a 1-chain, representing a line, into the two endpoints of that line. This operator has nilpotency, meaning the boundary operator applied to that same geometry is always zero. This property is important in the exterior derivative in preserving geometric invariants and physical conservation laws. 

This operator is computed by comparing the orientation of the k-chain with its respective (k-1)-chain boundary. It is represented as a matrix with dimensions $|C_(k−1) |×|C_k |$, meaning it inputs a k-chain and outputs a (k-1)-chain, as explained before. This matrix is sparse with elements of either 0, +1, or -1, depending on boundary relations with their respective k-chain. This computation only needs local information involving the k-chain and the simplexes within/around it. It does not require global information of the entire mesh to compute it. The computation of this operator is local, topological, and coordiante-free. As the geometric intuition for this matrix calculation is difficult to describe, I recommend further investigation online.

This tool is important in defining the exterior derivative operator through Generalized Stokes' Theorem. This theorem states that the exterior derivative of a form integrated/evaluated over a manifold is equivalent to the form integrated over the boundary of that manifold. Using this definition and the boundary operator's duality, the exterior derivative is the transpose of the boundary matrix. 

Finally, this operator is natural with respect to discrete pullbacks. This statement means that applying the boundary to simplex X before mapping a form on X to Y is the same as mapping the form on X to Y and then applying the operator to Y. This trait preserves the structure of the mesh under situations like mesh refinement and mesh deformations.

### Important Properties

∂(∂c)=0 

The formula above displays the property of nilpotency for the boundary operator. In other words, the boundary of a boundary is always zero. For example, the boundary of a line is its two endpoints, and the boundary of points is zero. The boundary of a ball is a hollow sphere, and the boundary of a hollow sphere is zero. This property applies to discrete chains as well.

⟨dω, c⟩ = ⟨ω, ∂c⟩

The boundary operator is dual to the exterior derivative. In fact, the exterior derivative can be called the coboundary operator. In other words, the exterior derivative of a form dω applied to the chain c is equivalent to the regular form ω applied to the boundary of that chain ∂c. This property is described in Generalized Stokes' Theorem and helps in the calculation of the exterior derivative. Consider ⟨ ⟩  to be brackets representing the inner product, where the inner product of a chain and a cochain is synonymous to integrating the cochain on the chain, requiring neither a metric nor a dot product. It is also the discrete analogue of integrating a form over a region.

$\partial[v_o, v_1, ..., v_k]=\sum^k_{i=0}(-1)^i[v_o, ..., v_i, ..., v_k]$

The formula above displays how the boundary operator is computed/defined in the discrete setting. The operator inputs a chain with a particular ordering of vertices, where edges may look like $e=[v_0, v_1]$ and faces may look like $σ=[v_0, v_1, v_2]$. Let's assume the operator inputs the triangular face $σ_o=[v_0, v_1, v_2]$ surrounded by the edges $e_0=[v_0, v_1]$, $e_1=[v_1, v_2]$, and $e_2=[v_o, v_2]$. The boundary operator omits the i-th vertex of the face for each iteration. The operation determines the sign of this edge based on if 'i' is even or odd. For this face $σ_o$, the result $∂[v_o, v_1, v_2 ]=[v_1, v_2 ]−[v_o, v_2 ]+[v_o, v_1 ]=e_1−e_2+e_o$. If the reader creates a visual representation of the face and surrounding vertices, with orientations going from the first vertex and moving to the last, the result will be the same.

### Citations

Discrete Exterior Calculus - Hirani  
Section 3.6 (pp. 35-36)

Discrete Exterior Calculus - Desbrun  
Section 5 (pp. 13-14)

Notes on Discrete Exterior Calculus  
Section 2.6 (p. 10)

Discrete Differential Forms  
Section 3 

## Exterior Derivative

### Signature

Decapodes.jl
    
    :d₀ => dec_differential(0, sd) |> matmul
    :d₁ => dec_differential(1, sd) |> matmul
    :dual_d₀ || :d̃₀ => dec_dual_derivative(0, sd) |> matmul
    :dual_d₁ || :d̃₁ => dec_dual_derivative(1, sd) |> matmul

    # Dual Differentials
    :dual_d₀ || :d̃₀ => dec_dual_derivative(0, sd) |> matmul
    :dual_d₁ || :d̃₁ => dec_dual_derivative(1, sd) |> matmul

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    d(s::HasDeltaSet, x::DualForm{n}) where n
    dual_derivative(n::Int, s::HasDeltaSet, args...) dual_derivative(::Type{Val{n}}, s::HasDeltaSet, args...) where n
    dual_derivative_nz(::Type{Val{0}}, s::AbstractDeltaDualComplex1D, x::Int) dual_derivative_nz(::Type{Val{0}}, s::AbstractDeltaDualComplex2D, x::Int) dual_derivative_nz(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, x::Int)

CombinatorialSpaces.jl, FastDEC.jl

    dec_differential(n::Int, sd::HasDeltaSet)
    dec_dual_derivative(n::Int, sd::HasDeltaSet)
    dec_p_dual_derivative(::Type{Val{0}}, sd::HasDeltaSet1D)
    dec_p_dual_derivative(::Type{Val{0}}, sd::HasDeltaSet2D)
    dec_p_dual_derivative(::Type{Val{1}}, sd::HasDeltaSet2D)
    dec_p_dual_derivative(::Type{Val{0}}, sd::HasDeltaSet3D)
    dec_p_dual_derivative(::Type{Val{1}}, sd::HasDeltaSet3D)
    dec_p_dual_derivative(::Type{Val{2}}, sd::HasDeltaSet3D)
    dec_p_derivbound(::Type{Val{0}}, sd::HasDeltaSet; transpose::Bool=false, negate::Bool=false)
    dec_p_derivbound(::Type{Val{1}}, sd::HasDeltaSet; transpose::Bool=false, negate::Bool=false)
    dec_p_derivbound(::Type{Val{2}}, sd::HasDeltaSet; transpose::Bool=false, negate::Bool=false)

### Description

The exterior derivative operator (ED) extends the idea of differentiation into differential geometry. It generalizes vector calculus concepts such as gradient, divergence, and curl onto multidimensional forms. It is a requirement for modelling any partial differential equation in an exterior calculus framework.

This operator maps k-forms to (k+1)-forms (or k-cochains to (k+1)-cochains in the discrete setting), representing the form's rate of change on the mesh. In 3D space, the ED operator represents the gradient on a 0-form, the curl on a 1-form, and divergence on a 2-form. The ED on a 3-form in 3D space is always zero.

This operator matrix is calculated with Generalized Stokes' Theorem, which states that the exterior derivative of a form integrated/evaluated over a manifold is equivalent to the form integrated over the boundary of that manifold. In the discrete setting, the exterior derivative for a cochain is the transpose of the boundary matrix operator for its respective chain. Like the boundary operator, it is a sparse matrix with elements 0, +1, and -1. The computation of this operator is local, topological, and coordiante-free.

Using geometric intuition, this operator determines the dimensional rate of change of a form based on the values of the forms around it. For example, imagine three cochains attatched to 1-simplexes (lines) that form a triangle. The exterior derivative of the 1-cochain formed from the three values is a 2-cochain residing on the triangular face, calculated by adding the cochains together and accounting for differing orientation.

Furthermore, the exterior derivative shares the property of nilpotency with its dual operator. Thus, the exterior derivative applied twice to the same form is always zero. This trait ingrains geometrical information and mathematical conservation properties, such as the curl of a form's divergence always being zero. 

On a dual mesh, the discrete exterior derivative may change in sign to account for the mesh's differing orientation. This formula is displayed in the Important Properties section.

Finally, this operator is natural with respect to discrete pullbacks, like the boundary. It also has an adjoint with the codifferential operator, allowing for Laplacians and Hodge decompositions.

### Important Properties

k-form $-->$ (k+1)-form

$d(dω)=0$

This is the property that the exterior derivative applied twice to the same form is always zero. This property retains the mathematical properties of $∇×∇f=0$ and $∇ ⋅∇×f=0$ without requiring any additional features. 

$d(ω∧β)=dω∧β+(−1)^{deg⁡(ω)}  ω∧β$

The formula above displays the Leibniz Principle of exterior derivatives. It shows how the operator interacts with the wedge product.

$d=∂^T$

In the discrete setting, the exterior derivative matrix is computed as the transpose of the boundary operator. It inputs k-forms and outputs the exterior derivative values of (k+1)-forms.

$d_{n−k}^{Dual}=(−1)^k (d_{k−1}^{Primal} )^T$

Due to the orientation changes in dual meshes, the exterior derivative applied to dual forms must be altered to account for differing orientation and orthogonality. The formula above displays this correction.

$\int_\Omega d\omega=\int_{\partial\Omega}\omega$

The formula above is Generalized Stokes' Theorem. It states that a person can determine the exterior derivative of a form across a manifold by determining the value of that across the boundary of the manifold. As the exterior derivative generalizes the gradient, divergence, and the curl. Generalized Stokes' Theorem combines Stokes' Theorem, Green's Theorem, and the Divergence Theorem into one expression.

### Citations

Discrete Exterior Calculus - Hirani  
Section 3.6 (pp. 35-37)

Discrete Exterior Calculus - Desbrun  
Section 5 (pp. 13)

Notes on Discrete Exterior Calculus  
Section 2.6 (p. 10-11)

Discrete Differential Forms  
Section 4.1-4.3

## Wedge Product

### Signature

Decapodes.jl

    :∧₀₁ => dec_pair_wedge_product(Tuple{0,1}, sd)  # Returns (in-place kernel, out-of-place)
    :∧₁₀ => dec_pair_wedge_product(Tuple{1,0}, sd)
    :∧₀₂ => dec_pair_wedge_product(Tuple{0,2}, sd)
    :∧₂₀ => dec_pair_wedge_product(Tuple{2,0}, sd)
    :∧₁₁ => dec_pair_wedge_product(Tuple{1,1}, sd)
    
    :∧ᵖᵈ₁₁ => dec_wedge_product_pd(Tuple{1,1}, sd)
    :∧ᵖᵈ₀₁ => dec_wedge_product_pd(Tuple{0,1}, sd)
    :∧ᵈᵖ₁₁ => dec_wedge_product_dp(Tuple{1,1}, sd)
    :∧ᵈᵖ₁₀ => dec_wedge_product_dp(Tuple{1,0}, sd)
    
    :∧ᵈᵈ₁₁ => dec_wedge_product_dd(Tuple{1,1}, sd)
    :∧ᵈᵈ₁₀ => dec_wedge_product_dd(Tuple{1,0}, sd)
    :∧ᵈᵈ₀₁ => dec_wedge_product_dd(Tuple{0,1}, sd)

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ∧(s::HasDeltaSet, α::SimplexForm{k}, β::SimplexForm{l}) where {k,l}
    ∧(k::Int, l::Int, s::HasDeltaSet, args...)
    ∧(::Type{Tuple{k,l}}, s::HasDeltaSet, α, β) where {k,l}
    ∧(::Type{Tuple{0,0}}, s::HasDeltaSet, f, g, x::Int)
    ∧(::Type{Tuple{k,0}}, s::HasDeltaSet, α, g, x::Int) where k
    ∧(::Type{Tuple{0,k}}, s::HasDeltaSet, f, β, x::Int) where k
    ∧(::Type{Tuple{1,1}}, s::HasDeltaSet2D, α, β, x::Int)
    ∧(::Type{Tuple{2,1}}, s::HasDeltaSet3D, α, β, x::Int)
    ∧(::Type{Tuple{1,2}}, s::HasDeltaSet3D, α, β, x::Int)
    wedge_product_zero(::Type{Val{k}}, s::HasDeltaSet, f, α, x::Int) where k

CombinatorialSpaces.jl, FastDEC.jl

    ∧(s::HasDeltaSet, α::SimplexForm{k}, β::SimplexForm{l}) where {k,l}
    ∧(s::HasDeltaSet, α::SimplexForm{1}, β::DualForm{1})
    ∧(s::HasDeltaSet, α::DualForm{1}, β::SimplexForm{1})
    dec_wedge_product(::Type{Tuple{m,n}}, sd::HasDeltaSet, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing) where {m,n}
    dec_wedge_product(m::Int, n::Int, sd::HasDeltaSet)
    dec_wedge_product(::Type{Tuple{0,0}}, sd::HasDeltaSet, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing)
    dec_wedge_product(::Type{Tuple{k,0}}, sd::HasDeltaSet, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing) where {k}
    dec_wedge_product(::Type{Tuple{0,k}}, sd::HasDeltaSet, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing) where {k}
    dec_wedge_product(::Type{Tuple{1,1}}, sd::HasDeltaSet2D, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing)
    dec_wedge_product(::Type{Tuple{1,2}}, sd::HasDeltaSet3D, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing)
    dec_wedge_product(::Type{Tuple{2,1}}, sd::HasDeltaSet3D, backend=Val{:CPU}, arr_cons=identity, cast_float=nothing)
    dec_wedge_product_dd(::Type{Tuple{m,n}}, sd::HasDeltaSet) where {m,n}
    dec_wedge_product_dd(::Type{Tuple{0,1}}, sd::HasDeltaSet)
    dec_wedge_product_dd(::Type{Tuple{1,0}}, sd::HasDeltaSet)
    dec_wedge_product_dp(::Type{Tuple{m,n}}, sd::HasDeltaSet) where {m,n}
    dec_wedge_product_dp(::Type{Tuple{1,0}}, sd::HasDeltaSet)
    dec_wedge_product_dp(::Type{Tuple{1,1}}, sd::HasDeltaSet)
    dec_wedge_product_pd(::Type{Tuple{m,n}}, sd::HasDeltaSet) where {m,n}
    dec_wedge_product_pd(::Type{Tuple{0,1}}, sd::HasDeltaSet)
    dec_wedge_product_pd(::Type{Tuple{1,1}}, sd::HasDeltaSet)
    dec_c_wedge_product!(::Type{Tuple{j,k}}, res, α, β, p, c) where {j,k}
    dec_c_wedge_product(::Type{Tuple{m,n}}, α, β, wedge_cache) where {m,n}

### Description

The wedge product operator is used to "combine" a k-form and an l-form into a higher degree (k+l)-form. As an example, imagine two vector-like 1-cochains. The wedge product would be the signed area formed by dragging one of the forms across the other, like a determinant. For this reason, the wedge product has nilpotency, where the wedge product of the same form $(ω∧ω=0)$ is likely zero.

It is an anticommutative operator, meaning $β∧ω$ is often equal to $−ω∧β$ under most circumstances. It is also associative, meaning $(ω∧β)∧γ=ω∧(β∧γ)$. In the discrete case however, the associative property is limitted to only closed forms. The exterior derivative of closed forms is always zero. Later operators which use this tool as a component also face this limitation.

Another limitation is the wedge product's reliance on the metric. Certain operator formulas may require defined angles for calculating the result of this operator. However, other formulas may lack this limitation. Formulas that are not metric-dependent are natural with respect to discrete pullbacks.

Finally, the Leibniz Rule displays an important relationship between the wedge product and exterior derivative. This rule is approximated using the operator formulas. 

### Important Properties

$ deg⁡(ω)+deg⁡(β)=deg⁡(ω∧β) $

The wedge product combines a k-form and a p-form and creates a (k+p)-form.


$ 
\langle \omega^k \wedge \beta^l, \sigma^{k+l}\rangle := 
\frac{1}{(k+l)!} 
\sum_{r\in S_{k+l+1}} sign(\tau) \frac{|\sigma^{k+l} \cap \star v_{\tau(k)}|}{|\sigma^{k+l}|} 
\langle \omega, [v_{\tau(0)}, ..., v_{\tau(k)}] \rangle 
\langle \beta, [v_{\tau(k)}, ..., v_{\tau(k+l)}] \rangle 
$

This is the metric formulation for the discrete primal-primal wedge product operator, where ω is a is a discrete primal k-form and β is a discrete primal l-form. The LHS states that the evaluation of the wedge product of two forms on the appropraite face is equal to the RHS. The RHS states that for each permutation of the simplex $σ^{k+l}$  (where there are $(k+l+1)!$  orderings of the vertices), it finds the sign of that permutation $(sign(τ))$ and multiplies it by a geometric weight $(|σ^{k+l}∩⋆v_τ (k)| / |σ^(k+l) | )$. The weight determines the volume of the new simplex that is shared by the simplices touching the k-th vertex of the current permutation, divided by the volume of the new simplex. Finally, it multiplies this weight by the inner product of ω on a k-simplex formed from initial vertices $(⟨ω, [v_τ(0) , …, v_τ(k) ]⟩)$ and the inner product of β on a l-simplex formed from the vertices afterwards $(⟨β, [v_τ(k) , …, v_τ(k+l) ]⟩)$. The contributions from each permutation are summed together and normalized based on the forms' dimensions $(1/(k+l)!)$. As this formula relies on simplex magnitudes, it is metric.

$ 
\langle \omega^k \wedge \beta^l, \sigma^{k+l}\rangle := 
\frac{1}{(k+l+1)!} 
\sum_{r\in S_{k+l+1}} sign(\tau)  
\langle \omega, [v_{\tau(0)}, ..., v_{\tau(k)}] \rangle 
\langle \beta, [v_{\tau(k)}, ..., v_{\tau(k+l)}] \rangle 
$

The formula above is the topological primal-primal wedge product operator. It follows the same trend as the metric case. However, this one does not use a geometric weighing factor, and instead accounts for weight with normalization $(1/(k+l+1)!)$. 

$ ω∧β=(−1)^{deg⁡(ω)deg⁡(β)}⁡(β∧ω) $

This formula displays the anti-commutative property of wedge product. Switching the order of the forms on the wedge product may swap the sign of the resulting higher degree form. This is the extension of the idea that forms are antisymmetric tensors, and swapping the forms that create a higher form change the orientation of the area/volume/length.

$ (ω∧β)∧γ=ω∧(β∧γ) $

This formula displays the associative property of the wedge product operator in the smooth case. In the discrete case, this is only true for closed forms. The exterior derivative of a closed form is zero.

$ ω∧ω=0 $

This formula shows the nilpotency property of wedge product in the smooth case. Consider two differing 1-forms, ω and β, which can be considered dual to vectors. Then, think of the wedge product as an operator which finds the signed area of the parallelogram created by smearing one of those vectors across the edge of the other. If those two forms were in the same direction, like two ω forms, the smearing would only create a line with an area of zero.

### Citations

Discrete Exterior Calculus - Hirani  
Section 7.1-7.3 (pp. 71-76)

Discrete Exterior Calculus - Desbrun  
Section 8 (pp. 17-21)

## Hodge Star

### Signature

Decapodes.jl

    :⋆₀ => dec_hodge_star(0, sd, hodge=hodge) |> matmul
    :⋆₁ => dec_hodge_star(1, sd, hodge=hodge) |> matmul
    :⋆₂ => dec_hodge_star(2, sd, hodge=hodge) |> matmul
    
    :⋆₀⁻¹ => dec_inv_hodge_star(0, sd, hodge) |> matmul
    :⋆₁⁻¹ => dec_pair_inv_hodge(Val{1}, sd, hodge)  # Special: returns (in-place, out-of-place) solver
    :⋆₂⁻¹ => dec_inv_hodge_star(1, sd, hodge) |> matmul

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ⋆(s::HasDeltaSet, x::SimplexForm{n}; kw...) where n
    ⋆(n::Int, s::HasDeltaSet, args...; kw...)
    ⋆(::Type{Val{n}}, s::HasDeltaSet; hodge::DiscreteHodge=GeometricHodge()) where n
    ⋆(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector; hodge::DiscreteHodge=GeometricHodge()) where n
    ⋆(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector, ::DiagonalHodge) where n
    ⋆(::Type{Val{n}}, s::HasDeltaSet, ::DiagonalHodge) where n
    ⋆(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, ::GeometricHodge)
    ⋆(::Type{Val{0}}, s::AbstractDeltaDualComplex2D, ::GeometricHodge)
    ⋆(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, ::GeometricHodge)
    ⋆(::Type{Val{0}}, s::AbstractDeltaDualComplex2D, form::AbstractVector, ::GeometricHodge)
    ⋆(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, form::AbstractVector, ::GeometricHodge)
    ⋆(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, form::AbstractVector, ::GeometricHodge)
    ⋆(::Type{Val{n}}, s::AbstractDeltaDualComplex1D, ::GeometricHodge) where n
    ⋆(::Type{Val{n}}, s::AbstractDeltaDualComplex1D, form::AbstractVector, ::GeometricHodge) where n
    inv_hodge_star(n::Int, s::HasDeltaSet, args...; kw...)
    inv_hodge_star(::Type{Val{n}}, s::HasDeltaSet; hodge::DiscreteHodge=GeometricHodge()) where n
    inv_hodge_star(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector; hodge::DiscreteHodge=GeometricHodge()) where n
    inv_hodge_star(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector, ::DiagonalHodge) where n
    inv_hodge_star(::Type{Val{n}}, s::HasDeltaSet, ::DiagonalHodge) where n
    inv_hodge_star(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, ::GeometricHodge)
    inv_hodge_star(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, form::AbstractVector, ::GeometricHodge)
    inv_hodge_star(::Type{Val{0}}, s::AbstractDeltaDualComplex2D, ::GeometricHodge)
    inv_hodge_star(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, ::GeometricHodge)
    inv_hodge_star(::Type{Val{0}}, s::AbstractDeltaDualComplex2D, form::AbstractVector, ::GeometricHodge)
    inv_hodge_star(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, form::AbstractVector, ::GeometricHodge)
    inv_hodge_star(::Type{Val{n}}, s::AbstractDeltaDualComplex1D, ::GeometricHodge) where n
    inv_hodge_star(::Type{Val{n}}, s::AbstractDeltaDualComplex1D, form::AbstractVector, ::GeometricHodge) where n
    
CombinatorialSpaces.jl, FastDEC.jl

    dec_hodge_star(n::Int, sd::HasDeltaSet; hodge=GeometricHodge())
    dec_hodge_star(n::Int, sd::HasDeltaSet, ::DiagonalHodge)
    dec_hodge_star(n::Int, sd::HasDeltaSet, ::GeometricHodge)
    dec_hodge_star(::Type{Val{k}}, sd::HasDeltaSet, ::DiagonalHodge) where {k}
    dec_hodge_star(::Type{Val{j}}, sd::EmbeddedDeltaDualComplex1D, ::GeometricHodge) where {j}
    dec_hodge_star(::Type{Val{j}}, sd::EmbeddedDeltaDualComplex2D, ::GeometricHodge) where {j}
    dec_hodge_star(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, point_type}, ::GeometricHodge) where {float_type, point_type}
    dec_p_hodge_diag(::Type{Val{0}}, sd::EmbeddedDeltaDualComplex1D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex1D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{0}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{2}}, sd::EmbeddedDeltaDualComplex2D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{0}}, sd::EmbeddedDeltaDualComplex3D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex3D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{2}}, sd::EmbeddedDeltaDualComplex3D{Bool, float_type, _p} where _p) where float_type
    dec_p_hodge_diag(::Type{Val{3}}, sd::EmbeddedDeltaDualComplex3D{Bool, float_type, _p} where _p) where float_type
    dec_inv_hodge_star(n::Int, sd::HasDeltaSet; hodge=GeometricHodge())
    dec_inv_hodge_star(n::Int, sd::HasDeltaSet, ::DiagonalHodge)
    dec_inv_hodge_star(n::Int, sd::HasDeltaSet, ::GeometricHodge)
    dec_inv_hodge_star(::Type{Val{k}}, sd::HasDeltaSet, ::DiagonalHodge) where {k}
    dec_inv_hodge_star(::Type{Val{j}}, sd::EmbeddedDeltaDualComplex1D, ::GeometricHodge) where {j}
    dec_inv_hodge_star(::Type{Val{j}}, sd::EmbeddedDeltaDualComplex2D, ::GeometricHodge) where {j}
    dec_inv_hodge_star(::Type{Val{1}}, sd::EmbeddedDeltaDualComplex2D, ::GeometricHodge)
    dec_inv_hodge_star(::Type{Val{0}}, sd::EmbeddedDeltaDualComplex3D, ::GeometricHodge)
    dec_inv_hodge_star(::Type{Val{3}}, sd::EmbeddedDeltaDualComplex3D, ::GeometricHodge)
dec_inv_hodge_star(::Type{Val{j}}, sd::EmbeddedDeltaDualComplex3D, ::GeometricHodge) where {j}

### Description

The hodge star operator maps k-forms to their respective "orthogonal" (n-k)-forms based on the Riemannian metric (or inner product). 'n' is the number of dimensions in the space. It relates the primal mesh to the dual mesh and is an important component in later operators. For example, in 2D space, points are outputted as areas, lines become orthogonal lines, and areas become midpoints of the space. 

The discrete hodge star is computed using a formula which states that the original form evaluated on the original simplex divided by the simplex's magnitude is equivalent to the starred form evaluated on the "orthogonal" simplex divided by the "orthogonal" simplex's magnitude. This formula ensures dual and primal meshes contain the same physical information. The operator is also an isomorphism, meaning the hodge star applied twice to a form will return the original form (with a potential sign change depending on dimensions). 

The operator matrix is calculated through either two methods - the geometric hodge or the diagonal hodge. The geometric hodge is calculated as a sparse, diagonal matrix with elements dependent on the ratio between simplex and orthogonal simplex magnitudes. The matrix dimensions are $|C_k |×|C_k |$. The diagonal hodge is calculated faster but is less accurate.

This operator is metric-dependent. It is reliant on a metric to determine angles and orthogonality. The result of this operator may differ even if simplex connectivity, cochain values, and orientations remain the same. For example, depending on if the dual mesh is barycentric-centered or circumcentric-centered, the hodge star may return differing results. Furthermore, stretching the mesh would cause changes in simplex magnitudes on primal and dual meshes, thus affecting the output values from the hodge star operator.

Finally, as this operator is metric-dependent, it does not commute fully with respect to natural pullbacks. However, the code attempts to approximate this property as much as possible for accurate simulations.

### Important Properties

k-form → (n-k)-form

$⟨ω^k, β^k ⟩v=ω^k∧⋆β^k$

The formula above displays the identity of the hodge star in the smooth case. The inner product of  two k-dimensional forms multiplied by their volume element is equal to the wedge product between the first form and the hodge star of the second form. The volume element encodes the metric and orientation, ensuring that the LHS has the same dimensions and sign as the righthand side. Due to this definition, it is inherantly a metric-dependent operator.

$ 
\frac{1}{|⋆σ^k |}  
⟨\starω, ⋆σ^k ⟩ ≔ 
\frac{1}{|σ^p |}  ⟨ω, σ^k ⟩
$

The formula above is the definition of the discrete hodge star for $1≤k≤n−1$ forms. It states that the hodge star of a form integrated on the dual of a simplex and divided by the length of that dual simplex must be equal to the form evaluated on the primal simplex, divided by the length of that primal simplex. The dual form maintains the same physical information while transitioning to orthogonality. Forms that do not conform to the inequality have an additional term that changes sign to account for simplex orientation.

$⋆⋆ω=(−1)^k(n−k)ω $

The hodge star operator is an isomorphism. Using the hodge star on a form twice returns it to its original state, with a change in sign depending on dual orientation.

### Citations

Discrete Exterior Calculus - Hirani  
Section 4.1 (pp. 40-42)

Discrete Exterior Calculus - Desbrun  
Section 6 (pp. 14-15)

Notes on Discrete Exterior Calculus  
Section 2.9 (pp. 12-13)

Discrete Differential Forms  
Section 5.3-5.4

## Flat

### Signature

### Description

### Important Properties

### Citations


## Sharp

### Signature

### Description

### Important Properties

### Citations



## Codifferential

### Signature

### Description

### Important Properties

### Citations


## Interior Product

### Signature

### Description

### Important Properties

### Citations


## Lie Derivative

### Signature

Decapodes.jl

    :L₀ => add_Lie_1D!(Val{0}, ...)  # d ∘ ι (interior then dual derivative)
    :L₁ => add_Lie_1D!(Val{1}, ...)  # d ∘ ι (same in 1D)
    :L₂ => add_Lie_2D!(Val{2}, ...)  # d ∘ ι (dual 2-form)
    
    :ℒ₁ => ℒ_dd(Tuple{1,1}, sd)      # Full dual-dual Lie: -(d∘ι + ι∘d)

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ℒ(s::HasDeltaSet, X♭::EForm, α::DualForm{n}; kw...) where n
    lie_derivative_flat(n::Int, s::HasDeltaSet, args...; kw...)
    lie_derivative_flat(::Type{Val{0}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    lie_derivative_flat(::Type{Val{1}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    lie_derivative_flat(::Type{Val{2}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)

CombinatorialSpaces.jl, FastDEC.jl

    ℒ_dd(::Type{Tuple{1,1}}, s::SimplicialSets.HasDeltaSet)

### Description

The Lie Derivative operator inputs a k-form and outputs a k-form representing the rate of change of the form as it moves in the direction of a vector field. It is analogous to the directional derivative in vector calculus. A common application for the lie derivative is the formula for advection. 

The Lie Derivative's algebraic definition is displayed in Cartan's Magic Formula, which uses exterior derivative and interior product operators. This smooth formula is used to define the discrete primal-dual lie derivative, where X is a discrete vector field on the primal mesh and ω is a dual k-form.

However, this method may not satisfy the Leibniz Property of Lie Derivatives, as it uses the discrete wedge product operator (within the interior product operator). This operator is only associative on closed forms. 

Hirani notes an alternative flow-out definition in Section 8.5 in his thesis, Discrete Exterior Calculus, which does not suffer from this lack of associativity. However, this definition is not covered here.

### Important Properties

k-form → k-form

$L_X ω=i_X (dω)+d(i_X ω)$

This formula displays Cartan's Magic Formula, or an algebraic definition of the lie derivative using the interior product and the exterior derivative operators. This formula helps define the discrete primal-dual lie derivative operator on a simplicial complex. X is a discrete primal vector field and ω is a dual p-form.


$L_X (ω∧β)=L_X (ω)∧β+ω∧L_X (β)$

This is the Leibniz Property of the lie derivative, which displays the lie derivative's relationship with wedge product in the smooth case. In the discrete case, if the algebraic definition of the interior product is used, this relationship might not always be correct. This is because the algebraic interior product is defined using the wedge product. As the wedge product is only associative on closed forms (in the discrete setting), there may be errors in this Leibniz property.

### Citations

Discrete Exterior Calculus - Hirani  
Section 8.4-8.5 (pp. 83-85)

Discrete Exterior Calculus - Desbrun  
Section 10 (pp. 25-26)

## Laplace-Beltrami

### Signature

Decapodes.jl

    :Δ₀ => add_De_Rham_1D!(Val{0}, ...)  # dδ + δd → expands to full Laplace–de Rham
    :Δ₁ => add_De_Rham_1D!(Val{1}, ...) or add_De_Rham_2D!(Val{1}, ...)
    :Δ₂ => add_De_Rham_2D!(Val{2}, ...)
    
    :Δᵈ₀ => Δᵈ(Val{0}, sd)               # Dual 0-form Laplacian (fast version)
    :Δᵈ₁ => Δᵈ(Val{1}, sd)               # Dual 1-form Laplacian
    
    :Δ₀⁻¹ => dec_inv_lap_solver(Val{0}, sd)  # factorize(∇²(0, sd)) \ x

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ∇²(s::HasDeltaSet, x::SimplexForm{n}; kw...) where n
    ∇²(n::Int, s::HasDeltaSet, args...; kw...)
    ∇²(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector; kw...) where n
    ∇²(::Type{Val{n}}, s::HasDeltaSet; matrix_type::Type=SparseMatrixCSC{Float64}, kw...) where n
    
    Δ(s::HasDeltaSet, x::SimplexForm{n}; kw...) where n
    Δ(n::Int, s::HasDeltaSet, args...; kw...)
    Δ(::Type{Val{0}}, s::HasDeltaSet, form::AbstractVector; kw...)
    Δ(::Type{Val{0}}, s::HasDeltaSet; matrix_type::Type=SparseMatrixCSC{Float64}, kw...)
    Δ(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector; kw...) where n
    Δ(::Type{Val{n}}, s::HasDeltaSet; matrix_type::Type=SparseMatrixCSC{Float64}, kw...) where n
    Δ(::Type{Val{1}}, s::AbstractDeltaDualComplex1D, form::AbstractVector; kw...)
    Δ(::Type{Val{1}}, s::AbstractDeltaDualComplex1D; matrix_type::Type=SparseMatrixCSC{Float64}, kw...)
    Δ(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, form::AbstractVector; kw...)
    Δ(::Type{Val{2}}, s::AbstractDeltaDualComplex2D; matrix_type::Type=SparseMatrixCSC{Float64}, kw...)
    
CombinatorialSpaces.jl, FastDEC.jl

    Δᵈ(::Type{Val{0}}, s::SimplicialSets.HasDeltaSet)
    Δᵈ(::Type{Val{0}}, s::SimplicialSets.HasDeltaSet2D)
    Δᵈ(::Type{Val{0}}, s::SimplicialSets.HasDeltaSet3D)
    Δᵈ(::Type{Val{1}}, s::SimplicialSets.HasDeltaSet)
    dec_Δ⁻¹(::Type{Val{0}}, s::AbstractGeometricMapSeries; scheme::AbstractSubdivisionScheme = BinarySubdivision(), steps = 3, cycles = 5, alg = cg, μ = 2)


### Description

The Laplace-Beltrami operator generalizes the traditional Laplacian to curved surfaces. It linearly maps k-forms to k-forms. The ordinary Laplacian determines how much a point in a function deviates from the average of the functions surrounding it. A positive Laplacian implies that the point is less than the surrounding average (or convex down), while a negative one implies the point is greater than the surrounding average (or convex up).

This operator is defined with the exterior derivative and codifferential operators. For 0-forms (scalars), this operator simplifies to $δd$, as the exterior codifferential of a 0-form is always zero.

It is represented as a sparse matrix. Due to its reliance on the codifferential (and thus the hodge star), it is metric-dependent. 

### Important Properties

k-form → k-form

$Δ=dδ+δd$

This is the general Laplace-deRham Operator, which converts k-forms to dimensionally equivalent k-forms in the smooth case. The Laplace-Beltrami operator, $∇f=δdf$, is a special situation for 0-form functions, where $d(δf)=0$ as the codifferential of a 0-form is always zero.

$⟨Δf, σ^0 ⟩ =\frac{1}{ |\sigma_o| }  ∑_{ σ^1=[σ^0, v] } \frac{|⋆σ^1 |}{|σ^1 |} (f(v) − f(σ^0)) $

This formula above showcases the evaluation of a 0-form (f) at a primal vertex (σ^0) on a well-orientated triangular mesh. The mesh does not need to be flat. The formula states that the Laplacian of a 0-form f at a vertex $σ^0$,  $(⟨Δf, σ^0 ⟩)$ can be evaluated by the sum/contribution of all edges bordering $σ^0$, each weighted by the length of their corresponding dual edge divded by their primal length $(|⋆σ^1 |/|σ^1 | )$, and then multiplied by the difference between the form at the evaluated vertex $(f(σ^0))$ and the neighboring vertex $(f(v))$, connected by the edge. Then, the result of the weighted sum is divided by the volume/area of the vertice's coresponding dual cell $(1/|⋆σ^0 | )$ to normalize it. This formula matches with another geometric Laplace-Beltrami formula that uses cotangents and indices.

### Citations

Discrete Exterior Calculus - Hirani  
Section 6.4 (pp. 68-70)

Discrete Exterior Calculus - Desbrun  
Section 10 (pp. 26-27)