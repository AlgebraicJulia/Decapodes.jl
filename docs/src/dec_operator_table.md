# Discrete Exterior Calculus Operator Table

This document is meant to provide simplified descriptions of the operators of Discrete Exterior Calculus. The reader can find more rigorous definitions in the referenced sources for each operator. The links to the sources referenced for the operators are below:

[Discrete Exterior Calculus](https://www.cs.jhu.edu/~misha/Fall09/Hirani03.pdf) - Hirani


[Discrete Exterior Calculus](https://arxiv.org/abs/math/0508341) - Desbrun et. al.

[Notes on Discrete Exterior Calculus](https://math.arizona.edu/~agillette/research/decNotes.pdf) - Gillette

[Discrete Differential Forms for Computational Modeling](https://geometry.caltech.edu/pubs/DKT05.pdf) - Desbrun et. al.  

## Boundary

### Signature

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

This operator inputs k-chains and outputs (k-1)-chains (or k-manifolds to (k-1)-manifolds in the smooth case). For example, it could map a 1-chain, representing a line, to the two endpoints of that line. This operator has nilpotency, meaning the boundary operator applied to the same geometry is always zero. This property is important in the exterior derivative for preserving geometric invariants and physical conservation laws. 

This operator is computed by comparing the orientation of the k-chain with its respective (k-1)-chain boundary. It is represented as a matrix with dimensions $|C_{k−1} |×|C_k |$. It inputs a k-chain and outputs a (k-1)-chain, as explained before. This matrix is sparse with elements of either 0, +1, or -1, depending on boundary relations with their respective k-chain. This computation only needs local information involving the k-chain and the simplexes within/around it. The computation of this operator is local, topological, and coordinate-free. As the geometric intuition for this matrix calculation is difficult to describe, I recommend further investigation online.

This tool is important in defining the exterior derivative operator through Generalized Stokes' Theorem. This theorem states that the exterior derivative of a form integrated/evaluated over a manifold is equivalent to the form integrated over the boundary of that manifold. Using this definition and the boundary operator's duality, the exterior derivative is the transpose of the boundary matrix. 

Finally, this operator is natural with respect to discrete pullbacks. This statement means that applying the boundary to simplex X before mapping a form on X to Y is the same as mapping the form on X to Y and then applying the operator to Y. This trait preserves the structure of the mesh under situations like mesh refinement and mesh deformations.

### Important Properties

$ \partial^k: Cₖ(K; ℤ) → Cₖ₋₁(K; ℤ) $  

The boundary operator linearly inputs linear combinations of k-chains (integers with respect to simplex orientation, represented by $ \Z $) on mesh K and outputs the combined (k-1)-chain sum of the boundary of those k-chains.  

$ ∂(∂c) = 0 $

The formula above displays the property of nilpotency for the boundary operator. In other words, the boundary of a boundary is always zero. For example, the boundary of a line is its two endpoints, and the boundary of points is zero. The boundary of a ball is a hollow sphere, and the boundary of a hollow sphere is zero. This property applies to chains as well.

$ ⟨dω, c⟩ = ⟨ω, ∂c⟩ $

The boundary operator is dual to the exterior derivative. In fact, the exterior derivative is defined as the coboundary operator. In other words, the exterior derivative of a form dω applied to the chain c is equivalent to the regular form ω applied to the boundary of that chain ∂c. This property is described in Generalized Stokes' Theorem and helps in the calculation of the exterior derivative. Consider ⟨ ⟩  to be brackets representing the inner product, where the inner product of a chain and a cochain is synonymous to integrating the cochain on the chain, not requiring a metric. It is also the discrete analogue of integrating a form over a region.

$ ∂[v₀, v₁, …, vₖ] = ∑ᵏᵢ₌₀ (−1)ⁱ [v₀, …, v̂ᵢ, …, vₖ] $  

The formula above displays how the boundary operator is computed/defined in the discrete setting. The operator inputs a chain with a particular ordering of vertices, where edges may look like $e=[v_0, v_1]$ and faces may look like $σ=[v_0, v_1, v_2]$. Let's assume the operator inputs the triangular face $σ_o=[v_0, v_1, v_2]$ surrounded by the edges $e_0=[v_0, v_1]$, $e_1=[v_1, v_2]$, and $e_2=[v_o, v_2]$. The boundary operator omits the i-th vertex of the face for each iteration. The operation determines the sign of this edge based on if 'i' is even or odd. For this face $σ_o$, the result $∂[v_o, v_1, v_2 ]=[v_1, v_2 ]−[v_o, v_2 ]+[v_o, v_1 ]=e_1−e_2+e_o$. If the reader creates a visual representation of the face and surrounding vertices, with orientations going from the first vertex and moving to the last, the result will be the same.

### Citations

Discrete Exterior Calculus - Hirani\
Section 3.6 (pp. 35-36)

Discrete Exterior Calculus - Desbrun et. al.\
Section 5 (pp. 13-14)

Notes on Discrete Exterior Calculus - Gillette\
Section 2.6 (p. 10)

Discrete Differential Forms for Computational Modeling - Desbrun et. al.\
Section 3 

## Exterior Derivative

### Signature

Decapodes.jl
    
    :d₀ => dec_differential(0, sd) |> matmul
    :d₁ => dec_differential(1, sd) |> matmul
    :dual_d₀ || :d̃₀ => dec_dual_derivative(0, sd) |> matmul
    :dual_d₁ || :d̃₁ => dec_dual_derivative(1, sd) |> matmul
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

This operator linearly maps k-forms to (k+1)-forms (or k-cochains to (k+1)-cochains in the discrete setting), representing the form's rate of change on the mesh. In 3D space, the ED operator represents the gradient on a 0-form, the curl on a 1-form, and divergence on a 2-form. The ED on a 3-form in 3D space is always zero.

This operator matrix is calculated with Generalized Stokes' Theorem, which states that the exterior derivative of a form integrated/evaluated over a manifold is equivalent to the form integrated over the boundary of that manifold. In the discrete setting, the exterior derivative for a cochain is the transpose of the boundary matrix operator for its respective chain. Like the boundary operator, it is a sparse matrix with elements 0, +1, and -1. The computation of this operator is local, topological, and coordinate-free.

This operator determines the dimensional rate of change of a form based on the values of the forms around it. For example, imagine three cochains attached to 1-simplexes (lines) that form a triangle. The exterior derivative of the 1-cochain formed from the three values is a 2-cochain residing on the triangular face, calculated by adding or subtracting the cochains together depending on face and edge orientations. 

Furthermore, the exterior derivative shares the property of nilpotency with its dual operator. Thus, the exterior derivative applied twice to the same form is always zero. This trait ingrains geometrical information and mathematical conservation properties, such as the curl of a form's divergence always being zero. 

On a dual mesh, the discrete exterior derivative may change in sign to account for the mesh's differing orientation. This formula is displayed in the Important Properties section.

Finally, this operator is natural with respect to discrete pullbacks, like the boundary. It also has an adjoint with the codifferential operator, allowing for Hodge Decompositions and the Laplacian.

### Important Properties

$ d^k: \Omega\_{d}^k (K) \rightarrow \Omega\_{d}^{k + 1} (K) $

The discrete exterior derivative linearly maps discrete k-forms on mesh K to discrete (k+1)-forms on mesh K.

$ d(dω)=0 $

The exterior derivative applied twice to the same form is always zero. This property retains the PDE's mathematical structure (geometric invariants, $∇×∇f=0$, and $∇ ⋅∇×f=0$) without requiring any additional features. 

$ d(ω∧β) = dω ∧ β + (−1)^{deg⁡(ω)} ω ∧ β $

The formula above displays the Leibniz Principle of exterior derivatives. The operator has a product-rule-like relationship with the wedge product.

$ d=∂^T $

In the discrete setting, the exterior derivative matrix is computed as the transpose of the boundary operator. It inputs k-forms and outputs the exterior derivative values of (k+1)-forms.

$ \tilde{d}\_{n−k} = (−1)^k (d\_{k−1} )^T $ 

Due to the orientation changes in dual meshes, the exterior derivative applied to dual forms must be altered to account for differing orientation and orthogonality. The formula above displays this correction, where $ \tilde{d} $ represents the dual exterior derivative operator.

$ \int \_{\Omega} d\omega = \int \_{\partial\Omega}\omega $ 

The formula above is Generalized Stokes' Theorem. It states that a person can determine the exterior derivative of a form across a manifold by determining the value of that across the boundary of the manifold. Thus, the exterior derivative generalizes the gradient, divergence, and the curl while Generalized Stokes' Theorem generalizes Stokes' Theorem, Green's Theorem, and the Divergence Theorem into one expression.

### Citations

Discrete Exterior Calculus - Hirani\
Section 3.6 (pp. 35-37)

Discrete Exterior Calculus - Desbrun et. al.\
Section 5 (pp. 13)

Notes on Discrete Exterior Calculus - Gillette\
Section 2.6 (p. 10-11)

Discrete Differential Forms for Computational Modeling - Desbrun et. al.\
Section 4.1-4.3

## Wedge Product

### Signature

Decapodes.jl

    :∧₀₁ => dec_pair_wedge_product(Tuple{0,1}, sd)  
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

The wedge product operator is used to "combine" a k-form and an l-form into a higher degree (k+l)-form. As an example, imagine two vector-like 1-cochains. The wedge product would be the signed area formed by dragging one of the forms across the other, like a determinant.

It is an anti-commutative operator, meaning $β∧ω$ is often equal to $−ω∧β$ under most circumstances. It is also associative, meaning $(ω∧β)∧γ=ω∧(β∧γ)$. In the discrete case however, the associative property is limited to only closed forms. The exterior derivative of closed forms is always zero. Later operators which use this tool as a component also face this limitation.

Another limitation is the wedge product's reliance on the metric. Certain operator formulas may require defined angles for calculating the result of this operator. However, other formulas may lack this limitation. Formulas that are not metric-dependent are natural with respect to discrete pullbacks.

Finally, the Leibniz Rule highlights an important relationship between the wedge product and exterior derivative. This rule is approximated using the operator formulas shown below. 

### Important Properties

$ \wedge:\Omega^k\_d (K) \times \Omega^l\_d (K) \rightarrow \Omega^{k + l}\_d (K) $

The discrete primal-primal wedge product inputs a discrete primal k-form and discrete primal l-form on mesh K and outputs a discrete (k+l)-form on mesh K. The discrete dual-dual wedge product and discrete primal-dual wedge product would be defined differently.

$ 
\langle \omega^k \wedge \beta^l, \sigma^{k+l} \rangle := 
\frac{1}{(k+l)!} 
\sum\_{r\in S\_{k+l+1}} sign(\tau) \frac{|\sigma^{k+l} \cap \star v\_{\tau(k)}|}{|\sigma^{k+l}|} 
\langle \omega, [v\_{\tau(0)}, ..., v\_{\tau(k)}] \rangle 
\langle \beta, [v\_{\tau(k)}, ..., v\_{\tau(k+l)}] \rangle 
$

This is the metric formulation for the discrete primal-primal wedge product operator, where ω is a is a discrete primal k-form and β is a discrete primal l-form. The LHS states that the evaluation of the wedge product of two forms on the appropriate face is equal to the RHS. The RHS states that for each permutation of the simplex $σ^{k+l}$  (where there are $(k+l+1)!$  orderings of the vertices), it finds the sign of that permutation $(sign(τ))$ and multiplies it by a geometric weight $(|σ^{k+l}∩⋆v_τ (k)| / |σ^(k+l) | )$. The weight determines the volume of the new simplex that is shared by the simplices touching the k-th vertex of the current permutation, divided by the volume of the new simplex. Finally, it multiplies this weight by the inner product of ω on a k-simplex formed from initial vertices $(⟨ω, [v_τ(0) , …, v_τ(k) ]⟩)$ and the inner product of β on a l-simplex formed from the vertices afterwards $(⟨β, [v_τ(k) , …, v_τ(k+l) ]⟩)$. The contributions from each permutation are summed together and normalized based on the forms' dimensions $(1/(k+l)!)$. As this formula relies on simplex magnitudes, it is metric.

$ 
\langle \omega^k \wedge \beta^l, \sigma^{k+l} \rangle := 
\frac{1}{(k+l+1)!} 
\sum\_{r\in S\_{k+l+1}} sign(\tau)  
\langle \omega, [v\_{\tau(0)}, ..., v\_{\tau(k)}] \rangle 
\langle \beta, [v\_{\tau(k)}, ..., v\_{\tau(k+l)}] \rangle 
$

The formula above is the topological primal-primal wedge product operator. It follows the same trend as the metric case. However, this one does not use a geometric weighing factor, and instead accounts for weight with normalization $(1/(k+l+1)!)$. 

$ ω∧β=(−1)^{deg⁡(ω)deg⁡(β)}⁡(β∧ω) $

This formula displays the anti-commutative property of wedge product. Switching the order of the forms on the wedge product may swap the sign of the resulting higher degree form. This is the extension of the idea that forms are antisymmetric tensors, and swapping the forms that create a higher form change the orientation of the area/volume/length.

$ (ω∧β)∧γ = ω∧(β∧γ) $

This formula displays the associative property of the wedge product operator in the smooth case. In the discrete case, this is only true for closed forms. The exterior derivative of a closed form is zero.

$ d(ω∧β) = dω ∧ β+(−1)^{deg⁡(ω)} ω ∧ β $

As displayed in the exterior derivative section, the Leibniz Rule defines the interaction between the exterior derivative and the wedge product operators. 

$ ω ∧ ω = 0 $

This formula shows the alternating property of wedge product in the smooth case. Due to its antisymmetric nature, the wedge product of two equivalent, even-dimensioned forms is always zero.

### Citations

Discrete Exterior Calculus - Hirani\
Section 7.1-7.3 (pp. 71-76)

Discrete Exterior Calculus - Desbrun et. al.\
Section 8 (pp. 17-21)

## Hodge Star

### Signature

Decapodes.jl

    :⋆₀ => dec_hodge_star(0, sd, hodge=hodge) |> matmul
    :⋆₁ => dec_hodge_star(1, sd, hodge=hodge) |> matmul
    :⋆₂ => dec_hodge_star(2, sd, hodge=hodge) |> matmul
    :⋆₀⁻¹ => dec_inv_hodge_star(0, sd, hodge) |> matmul
    :⋆₁⁻¹ => dec_pair_inv_hodge(Val{1}, sd, hodge)  
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

The hodge star operator maps k-forms to their respective "orthogonal" (n-k)-forms based on the Riemannian metric (or inner product). 'n' is the number of dimensions in the space. It relates the primal mesh to the dual mesh. For example, in 2D space, points are outputted as areas, lines become orthogonal lines, and areas become midpoints of the space. It is also an important component in later operators.

The discrete hodge star is computed using a formula which states that the original form evaluated on the original simplex divided by the simplex's magnitude is equivalent to the starred form evaluated on the "orthogonal" simplex divided by the "orthogonal" simplex's magnitude. This formula ensures dual and primal meshes contain the same information. The operator is also an isomorphism, meaning the hodge star applied twice to a form will return the original form (with a potential sign change depending on dimensions). 

The operator matrix is calculated through either two methods - the geometric hodge or the diagonal hodge. The geometric hodge is calculated as a sparse, symmetrical matrix with elements dependent on the ratio between simplex and orthogonal simplex magnitudes. It pairs well with barycentric dual mesh centering. The matrix dimensions are $|C_k |×|C_k |$. The diagonal hodge is a diagonal matrix with the same dimensions. It is calculated faster but has less accuracy. It pairs well with circumcentric dual mesh centering.

This operator is metric-dependent. It is reliant on a metric to determine angles and orthogonality. Thus, the result of this operator may differ even if simplex connectivity, cochain values, and orientations remain the same. For example, depending on if the dual mesh is barycentric-centered or circumcentric-centered, the hodge star may return differing results. Furthermore, stretching the mesh would cause changes in simplex magnitudes on primal and dual meshes, thus affecting the output values from the hodge star operator.

Finally, as this operator is metric-dependent, it does not commute fully with respect to natural pullbacks. However, the code attempts to approximate this property as much as possible for accurate simulations.

### Important Properties

$ *^k:\Omega^k\_d (K) \rightarrow \Omega^{n - k}\_d (\star K) $

The discrete hodge star linearly maps discrete k-forms on mesh K to discrete (n-k)-forms on dual mesh $ \star K$, where n is the number of dimensions of the space. The $ \star $ symbol acts on geometry and outputs its dual version.

$ ⟨ω^k, β^k ⟩v = ω^k∧*β^k $

The formula above displays the identity of the hodge star in the smooth case. The inner product of  two k-dimensional forms multiplied by their volume element is equal to the wedge product between the first form and the hodge star of the second form. The volume element encodes the metric and orientation, ensuring that the LHS has the same dimensions and sign as the right-hand side. Due to this definition, it is inherently a metric-dependent operator.

$ 
\frac{1}{|⋆σ^k |}  
⟨*ω, ⋆σ^k ⟩ ≔ 
\frac{1}{|σ^k |}  ⟨ω, σ^k ⟩
$

The formula above is the definition of the discrete hodge star for $1≤k≤n−1$ forms. It states that the hodge star of a form integrated on the dual of a simplex and divided by the length of that dual simplex must be equal to the form evaluated on the primal simplex, divided by the length of that primal simplex. The dual form maintains the same physical information while transitioning to orthogonality. Forms that do not conform to the inequality have an additional term that changes sign to account for simplex orientation.

$ **ω = (−1)^k (n−k)ω $

The hodge star operator is an isomorphism. Using the hodge star on a form twice returns it to its original state, with a change in sign depending on dual orientation.

### Citations

Discrete Exterior Calculus - Hirani\
Section 4.1 (pp. 40-42)

Discrete Exterior Calculus - Desbrun et. al.\
Section 6 (pp. 14-15)

Notes on Discrete Exterior Calculus - Gillette\
Section 2.9 (pp. 12-13)

Discrete Differential Forms for Computational Modeling - Desbrun et. al.\
Section 5.3-5.4

## Flat

### Signature

Decapodes.jl

    :♭ᵈᵖ => dec_♭(sd)                    
    :♭♯ => ♭♯_mat(sd) |> matmul          


CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ♭(s::HasDeltaSet, X::DualVectorField)
    ♭(s::AbstractDeltaDualComplex2D, X::AbstractVector, ::DPPFlat)
    ♭(s::AbstractDeltaDualComplex2D, X::AbstractVector, ::PPFlat)
    ♭_mat(s::AbstractDeltaDualComplex2D, f::DPPFlat)
    ♭_mat(s::AbstractDeltaDualComplex2D, p2s, ::DPPFlat)
    ♭_mat(s::AbstractDeltaDualComplex2D, ::PPFlat)


### Description

The flat operator is a musical isomorphism that linearly maps vector fields to 1-forms. In the smooth case, it is the inverse of the sharp operator. In the discrete case, it has no exact inverse. It is often used to define certain vector quantities as forms, such as fluid flow.

The smooth definition of the flat operator is reliant on an inner product, making the operator metric-dependent. Due to the existence of primal meshes, dual meshes, and differing interpolation methods, there are a total of eight different discrete flat operators. One example is the $♭_{dpp}$  operator, or the dual-primal-primal operator. This operator inputs discrete dual vector fields, uses dual-primal interpolation, and outputs cochains on the primal mesh. For non-flat meshes, operators that input dual vector fields are more applicable than ones that do not. This operator is represented as a sparse mesh.

### Important Properties

$ \flat: \mathfrak{X}\_d(K) \rightarrow \Omega^1\_d(K) $

The discrete flat operator maps discrete vector fields (usually tangent vectors at primal vertices or circumcenters/barycenters) on mesh K to discrete 1-forms on K.

$ ⟨X, v⟩ = X^♭ (v)$

The formula defines the flat operator in the smooth case. The inner product of a vector field X with a second vector field v is equivalent to the inner product $ X^♭ (v) = \langle \rangle $ at every point of the Riemannian manifold M. The metric is defined with $ \langle \rangle $.

$ \langle X^{\flat\_{dpp}}, \sigma^1 \rangle =
\sum_{\sigma^n \succ \sigma^1} \frac{| \star \sigma^1 \cap \sigma^n |}{| \star \sigma^1 |}
X(\sigma^n)\cdot \tilde{\sigma^1}
$

The formula above displays the definition for the DPP flat operator, or the operator that inputs a vector field on the dual mesh, uses dual-primal interpolation, and outputs a cochain on the primal mesh. The LHS states that the flatted vector field evaluated on an edge is equivalent to the RHS. The RHS states that for every simplex that contains the evaluated edge $(∑_{σ^n≻ σ^1} )$, determine how much of the evaluated edge's dual cell is within the simplex, find that magnitude, and then divide it by the magnitude of the edge's dual cell $( |\star σ^1 ∩ σ^n |/|\star σ^1 | )$. Then, multiply this value by the dot product of the simplex's average vector with the unit vector of the evaluated edge $(X(\sigma^n) \cdot \tilde{\sigma^1})$. Finally, sum all contributions.

### Citations

Discrete Exterior Calculus - Hirani\
Section 5.3-5.6 (pp. 46-54)

Discrete Exterior Calculus - Desbrun et. al.\
Section 7 (p. 16)

Notes on Discrete Exterior Calculus - Gillette\
Section 2.12 (pp. 14-15)

## Sharp

### Signature

Decapodes.jl

    :♯ᵖᵈ => dec_♯_pd(sd)                 
    :♯ᵖᵖ => dec_♯_pp(sd)                 
    :♯ᵈᵈ => dec_♯_dd(sd)                 

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ♯(s::HasDeltaSet2D, α::EForm)
    ♯(s::HasDeltaSet2D, α::DualForm{1})
    ♯(s::AbstractDeltaDualComplex1D, X::AbstractVector, ::PDSharp)
    ♯(s::AbstractDeltaDualComplex1D, X::AbstractVector, ::PPSharp)
    ♯(s::AbstractDeltaDualComplex2D, α::AbstractVector, DS::DiscreteSharp)
    ♯(s::AbstractDeltaDualComplex2D, α::AbstractVector, ::LLSDDSharp)
    ♯_mat(s::AbstractDeltaDualComplex2D, DS::DiscreteSharp)
    ♯_mat(s::AbstractDeltaDualComplex2D, ::LLSDDSharp)
    ♯_denominator(s::AbstractDeltaDualComplex2D, v::Int, t::Int, ::DiscreteSharp)
    ♯_denominator(s::AbstractDeltaDualComplex2D, v::Int, _::Int, ::AltPPSharp)
    get_orthogonal_vector(s::AbstractDeltaDualComplex2D, v::Int, e::Int)
    ♭♯(s::HasDeltaSet2D, α::SimplexForm{1})
    ♭♯_mat(s::HasDeltaSet2D)

### Description

The sharp operator is a musical isomorphism that linearly maps 1-forms to vector fields. In the smooth case, it is the inverse of the flat operator. In the discrete case, it has no exact inverse. It is used as a component in the interior product and Lie derivative operators.

As the smooth definition of the sharp is reliant on a metric, the sharp is a metric-dependent operator. There are a total of four different discrete sharp operators to account for form inputs and vector outputs on the primal and dual meshes. It is represented as a sparse matrix. 

The sharp vector turns k-forms into k-vectors. In the smooth case, it is inverses with the flat operator. The discrete sharp operator(s) lack this inverse.

### Important Properties

$ \sharp: \Omega^1\_d(K) \rightarrow \mathfrak{X}\_d(K) $

The discrete sharp operator maps discrete 1-forms on mesh K to discrete vector fields on mesh K. The vectors are usually tangent to primal vertices or mesh circumcenters/barycenters. 

$ ⟨ω^{\sharp}, v⟩ = ω(v) $

The formula defines the sharp operator map in the smooth case. The inner product $\omega^{\sharp}$ with vector field v is equivalent to the output of form $ \omega $ with the input v. This should be the case for all points in Riemannian manifold M.

### Citations

Discrete Exterior Calculus - Hirani\
Chapter 5.7-5.8 (pp. 54-56)

Discrete Exterior Calculus - Desbrun et. al.\
Section 7 (p. 16)

## Codifferential

### Signature

Decapodes.jl

    :δ₁ => add_Codiff!(d, src, tgt)  
    :δ₂ => add_Codiff!(d, src, tgt)

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    δ(s::HasDeltaSet, x::SimplexForm{n}; kw...) where n
    δ(n::Int, s::HasDeltaSet, args...; kw...)
    δ(::Type{Val{n}}, s::HasDeltaSet; hodge::DiscreteHodge=GeometricHodge(), matrix_type::Type=SparseMatrixCSC{Float64}) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector; hodge::DiscreteHodge=GeometricHodge()) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, ::DiagonalHodge, args...) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, ::GeometricHodge, matrix_type) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, ::GeometricHodge, form::AbstractVector) where n

### Description

The codifferential operator acts as the adjoint of the exterior derivative. It linearly maps k-forms to (k-1)-forms. This operator can represent divergence and it is used in the algebraic definition of the Laplacian operator.

This operator is metric due to the reliance of the hodge star in its definition. It is represented as a sparse matrix created by combining hodge star and exterior derivative matrices. Due to its simple definition, it can be computed just by applying three matrix operators and accounting for orientation.

Finally, due to the presence of the exterior derivative in its definition, the smooth codifferential operator has nilpotency - the operator applied to same form twice is zero. However, due to inaccuracies in the discrete hodge star, the discrete codifferential applied twice is only approximately zero.

### Important Properties

$ \delta^k: \Omega\_d^{k+1} (K) \rightarrow \Omega\_d^k (K) $

The discrete codifferential operator is a linear map that inputs discrete (k+1)-forms on mesh K and outputs discrete k-forms on K.

$ δf = 0 $

The formula above shows that the codifferential of a 0-form f is equal to zero. The codifferential operator takes a k-form and outputs a (k-1)-form. There is no such thing as a (-1)-form, so the operation outputs zero.

$ δω = (−1)^{n(k−1)+1} * d * ω $

The formula above displays the discrete definition of the codifferential operator on forms. It can be defined purely through hodge star and exterior derivative operators. The $(−1)^{n(k−1)+1}$  represents the orientation of the output. It is important to realize that if the codifferential is applied to a primal form, then the dual exterior derivative formula should be used in the definition, as $⋆ω$ would produce a dual form. 

$ ⟨dω, β⟩ = ⟨ω, δβ⟩ $

The exterior derivative and codifferential are adjoint in both the smooth and discrete settings. 

### Citations

Discrete Exterior Calculus - Hirani\
Section 4.2 (p. 42)

Discrete Exterior Calculus - Desbrun et. al.\
Section 6 (p. 15)

Discrete Differential Forms for Computational Modeling - Desbrun et. al.\
Section 5.5

## Interior Product

### Signature

Decapodes.jl

    :i₁ => add_Inter_Prod_1D! or add_Inter_Prod_2D!(Val{1}, ...)  
    :i₂ => add_Inter_Prod_2D!(Val{2}, ...)                         
    
    :ι₁₁ => interior_product_dd(Tuple{1,1}, sd)   
    :ι₁₂ => interior_product_dd(Tuple{1,2}, sd)   

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    interior_product(s::HasDeltaSet, X♭::EForm, α::DualForm{n}; kw...) where n
    interior_product_flat(n::Int, s::HasDeltaSet, args...; kw...)
    interior_product_flat(::Type{Val{n}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...) where n

CombinatorialSpaces.jl, FastDEC.jl

    interior_product_dd(::Type{Tuple{1,1}}, s::SimplicialSets.HasDeltaSet)
    interior_product_dd(::Type{Tuple{1,2}}, s::SimplicialSets.HasDeltaSet)


### Description

The interior product contracts a vector field with a k-form, outputting a (k-1)-form. It is an operator that "combines" forms and vector fields together by plugging in the field into one of the form's inputs. In exterior algebra, forms are measurement tools that input vectors and output scalar values. The interior product makes use of this definition to merge vector fields and forms together.

In the smooth setting, it is a topological operator. However, using the algebraic definition, this tool becomes metric due to its reliance on the flat and hodge star operators.

Due to the use of the wedge product in its definition, the Leibniz Property for this operator only applies for closed forms due to limits in associativity. This operator is used in the Lie derivative, which faces this limitation as well. 

Hirani notes the existence of an extrusion definition for this operator which is not covered in this description. However, it should be known that algebraic contraction is dual to the idea of extrusion in a geometric setting.

### Important Properties

$ i\_X: \Omega^k\_d (K) \rightarrow \Omega^{k-1}\_d (K) $

The discrete interior product operator $ i_X $ maps k-forms on discrete mesh K to (k-1)-forms on K, where X is a vector field.

$ i\_X ω(X\_1, …, X\_k ) = ω(X, X\_1. …, X\_k ) $

The above definition shows the algebraic definition of the interior product. This operator contracts a vector field with a k-form, "inserting" itself into one of the slots of the form. The output is a (k-1)-form.

$ i_X ω = (−1)^k(n−k) *(*ω∧X^♭ ) $

The above identity displays the definition of the interior product operator composed of hodge star, wedge product, and flat operators. The discrete operator can be created in the discrete setting with discrete forms of these operators. The wedge product, hodge star, and flat operations $(*ω∧X^{\flat})$ represents inserting the contribution of vector field X into the form ω. The exterior hodge star returns the result to its proper dimension, k+1. Finally, the $(−1)^k(n−k)$ accounts for changes in orientation due to the hodge star and wedge product operators. Due to the use of the discrete wedge product in this identity, the Leibniz Property for contraction will only apply to closed forms (where $dω=0$), as the discrete wedge product is only associative on closed forms.

$ i_X (f) = 0 $

The interior product of a 0-form f with a vector field X is always zero.

$ i\_X (ω^k ∧ β^l ) = (i\_X ω)∧β+(−1)^k ω∧(i\_X β) $

The interior product operator follows the Leibniz Rule. For non-closed forms, it is important to note that the discrete wedge product is not associative and thus errors in this rule may appear.

### Citations

Discrete Exterior Calculus - Hirani\
Section 8.2-8.3 (pp. 79-83)

Discrete Exterior Calculus - Desbrun et. al.\
Section 10 (pp. 24-26)

## Lie Derivative

### Signature

Decapodes.jl

    :L₀ => add_Lie_1D!(Val{0}, ...)  
    :L₁ => add_Lie_1D!(Val{1}, ...)  
    :L₂ => add_Lie_2D!(Val{2}, ...)  
    :ℒ₁ => ℒ_dd(Tuple{1,1}, sd)      

CombinatorialSpaces.jl, DiscreteExteriorCalculus.jl

    ℒ(s::HasDeltaSet, X♭::EForm, α::DualForm{n}; kw...) where n
    lie_derivative_flat(n::Int, s::HasDeltaSet, args...; kw...)
    lie_derivative_flat(::Type{Val{0}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    lie_derivative_flat(::Type{Val{1}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    lie_derivative_flat(::Type{Val{2}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)

CombinatorialSpaces.jl, FastDEC.jl

    ℒ_dd(::Type{Tuple{1,1}}, s::SimplicialSets.HasDeltaSet)

### Description

The Lie Derivative operator inputs a k-form and outputs a k-form representing the rate of change of the form as it moves in the direction of a vector field. It is analogous to the directional derivative in vector calculus. A common application for the Lie derivative is the formula for advection. 

The Lie derivative's algebraic definition is displayed in Cartan's Magic Formula, which uses exterior derivative and interior product operators. This smooth formula is used to define the discrete primal-dual lie derivative, where X is a discrete vector field on the primal mesh and ω is a dual k-form.

However, this method may not satisfy the Leibniz Property of Lie Derivatives, as it uses the discrete wedge product operator (within the interior product operator). This operator is only associative on closed forms. 

Hirani notes an alternative flow-out definition in Section 8.5 in his thesis, Discrete Exterior Calculus, which does not suffer from this lack of associativity. However, this definition is not covered here.

### Important Properties

$ \mathcal{L}\_X: \Omega^k\_d (K) \rightarrow \Omega^k\_d (K) $

The discrete Lie derivative $ \mathcal{L}_X $ maps k-forms on discrete mesh K to k-forms on K, where X is a vector field.

$ \mathcal{L}\_X ω = i\_X (dω) + d(i\_X ω) $

This formula displays Cartan's Magic Formula, or an algebraic definition of the Lie derivative using the interior product and the exterior derivative operators. This formula helps define the discrete primal-dual Lie derivative operator on a simplicial complex. X is a discrete primal vector field and ω is a dual p-form.


$ \mathcal{L}\_X (ω ∧ β) = \mathcal{L}\_X (ω) ∧ β + ω ∧ \mathcal{L}\_X (β) $

This is the Leibniz Property of the Lie derivative, which displays the Lie derivative's relationship with wedge product in the smooth case. In the discrete case, if the algebraic definition of the interior product is used, this relationship might not always be correct. This is because the algebraic interior product is defined using the wedge product. As the wedge product is only associative on closed forms (in the discrete setting), there may be errors in this Leibniz property.

### Citations

Discrete Exterior Calculus - Hirani\
Section 8.4-8.5 (pp. 83-85)

Discrete Exterior Calculus - Desbrun et. al.\
Section 10 (pp. 25-26)

## Laplacian Operator

### Signature

Decapodes.jl

    :Δ₀ => add_De_Rham_1D!(Val{0}, ...)  
    :Δ₁ => add_De_Rham_1D!(Val{1}, ...) or add_De_Rham_2D!(Val{1}, ...)
    :Δ₂ => add_De_Rham_2D!(Val{2}, ...)
    :Δᵈ₀ => Δᵈ(Val{0}, sd)               
    :Δᵈ₁ => Δᵈ(Val{1}, sd)              
    :Δ₀⁻¹ => dec_inv_lap_solver(Val{0}, sd)  

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

The Laplacian operator, or the Laplace-Beltrami operator on 0-forms, generalizes the traditional Laplacian onto curved surfaces. It linearly maps k-forms to k-forms. The ordinary Laplacian determines how much a point in a function deviates from the average of the points surrounding it. A positive Laplacian implies that the point is less than the surrounding average (or concave up), while a negative one implies the point is greater than the surrounding average (or concave up).

This operator is defined with the exterior derivative and codifferential operators. For 0-forms (scalars), this operator simplifies from $ \Delta = \delta d + d \delta $ to $ δd $, as the codifferential of a 0-form is always zero.

It is represented as a sparse matrix. Due to its reliance on the codifferential (and thus the hodge star), it is metric-dependent. 

### Important Properties

$ \Delta^k:\Omega^k\_d (K) \rightarrow \Omega^k\_d (K) $

The discrete Laplace-deRham operator, the general form of the Laplace-Beltrami, maps k-forms on discrete mesh K to k-forms on the same mesh.

$ Δ = dδ + δd $

This is the general Laplace-deRham Operator, which converts k-forms to dimensionally equivalent k-forms in the smooth case. The Laplace-Beltrami operator, $ \Delta f=δdf $, is a special situation for 0-form functions where $ d(δf) = 0$ as the codifferential of a 0-form is always zero.

$ ⟨Δf, σ^0 ⟩ = \frac{1}{ |\sigma\_o| }  ∑\_{ σ^1=[σ^0, v] } \frac{|⋆σ^1 |}{|σ^1 |} (f(v) − f(σ^0)) $

This formula above showcases the evaluation of a 0-form (f) at a primal vertex (σ^0) on a well-orientated triangular mesh. The mesh does not need to be flat. The formula states that the Laplacian of a 0-form f at a vertex $σ^0$,  $(⟨Δf, σ^0 ⟩)$ can be evaluated by the sum/contribution of all edges bordering $σ^0$, each weighted by the length of their corresponding dual edge divded by their primal length $(|⋆σ^1 |/|σ^1 | )$, and then multiplied by the difference between the form at the evaluated vertex $(f(σ^0))$ and the neighboring vertex $(f(v))$, connected by the edge. Then, the result of the weighted sum is divided by the volume/area of the vertice's coresponding dual cell $(1/|⋆σ^0 | )$ to normalize it. This formula was found to be equivalent with a different geometric Laplace-Beltrami formula that uses cotangents and indices.

### Citations

Discrete Exterior Calculus - Hirani\
Section 6.4 (pp. 68-70)

Discrete Exterior Calculus - Desbrun et. al.\
Section 10 (pp. 26-27)