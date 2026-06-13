# Discrete Exterior Calculus Operator Table

This document provides simplified descriptions of the operators of discrete exterior calculus (DEC). More rigorous definitions can be found in the cited sources for each operator. Links to those sources appear at the bottom of the document.

## Boundary

### Signature

#### CombinatorialSpaces.jl

    ∂(s::HasDeltaSet, x::DualChain{n}) where n
    dual_boundary(n::Int, s::HasDeltaSet, args...)
    dual_boundary(::Type{Val{n}}, s::HasDeltaSet, args...) where n
    dual_boundary_nz(::Type{Val{1}}, s::AbstractDeltaDualComplex1D, x::Int)
    dual_boundary_nz(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, x::Int)
    dual_boundary_nz(::Type{Val{2}}, s::AbstractDeltaDualComplex2D, x::Int)
    dec_boundary(n::Int, sd::HasDeltaSet)
    dec_p_boundary(::Type{Val{k}}, sd::HasDeltaSet; negate::Bool=false) where {k}

### Description

The boundary operator acts on the geometry of the mesh, linearly mapping ``k``-chains to ``(k-1)``-chains. It relates to the exterior derivative operator via duality and helps compute that operator through [Generalized Stokes' Theorem](https://en.wikipedia.org/wiki/Generalized_Stokes_theorem).

This operator inputs ``k``-chains and outputs ``(k-1)``-chains (or ``k``-manifolds to ``(k-1)``-manifolds in the smooth case). For example, it could map a 1-chain, representing a line, to the two endpoints of that line. This operator is nilpotent, meaning that applying the boundary operator twice to the same geometry always gives zero. This property is important in the exterior derivative for preserving geometric invariants and physical conservation laws.

This operator is computed by comparing the orientation of the ``k``-chain with its respective ``(k-1)``-chain boundary. It is represented as a matrix with dimensions ``|C_{k-1}| \times |C_k|``. It inputs a ``k``-chain and outputs a ``(k-1)``-chain, as explained before. This matrix is sparse with elements of either 0, +1, or -1, depending on boundary relations with their respective ``k``-chain. This computation only needs local information involving the ``k``-chain and the simplices incident to it. The computation of this operator is local, topological, and coordinate-free.

The boundary operator is important in defining the exterior derivative operator through Generalized Stokes' Theorem. This theorem states that integrating the exterior derivative of a form over a manifold is equivalent to integrating the form over the boundary of that manifold. Using this duality, the discrete exterior derivative is represented by the transpose of the boundary matrix.

Finally, this operator is natural with respect to discrete pullbacks. This means that applying the boundary to a simplex in ``X`` before applying a chain map from ``X`` to ``Y`` gives the same result as first applying the chain map and then taking the boundary of the mapped simplex. This property preserves the mesh structure under operations such as mesh refinement and mesh deformations.

```math
\begin{array}{ccc}
C_k(X) & \xrightarrow{f_\#} & C_k(Y) \\
\downarrow\scriptstyle{\partial_X} & & \downarrow\scriptstyle{\partial_Y} \\
C_{k-1}(X) & \xrightarrow{f_\#} & C_{k-1}(Y)
\end{array}
```

### Important Properties

```math
\partial^k: C_k(K; \mathbb{Z}) \to C_{k-1}(K; \mathbb{Z})
```

The boundary operator linearly inputs linear combinations of ``k``-chains (integers with respect to simplex orientation, represented by ``\mathbb{Z}``) on mesh ``K`` and outputs the combined ``(k-1)``-chain sum of the boundary of those ``k``-chains.

```math
\partial(\partial c) = 0
```

The formula above displays the property of nilpotency for the boundary operator. In other words, the boundary of a boundary is always zero. For example, the boundary of a line is its two endpoints, and the boundary of points is zero. The boundary of a ball is a hollow sphere, and the boundary of a hollow sphere is zero. This property applies to chains as well.

```math
\langle d\omega, c \rangle = \langle \omega, \partial c \rangle
```

The boundary operator is dual to the exterior derivative. In fact, the exterior derivative is defined as the coboundary operator. In other words, evaluating the exterior derivative ``d\omega`` on the chain ``c`` is equivalent to evaluating ``\omega`` on the boundary ``\partial c``. This relationship is captured by Generalized Stokes' Theorem and helps define the exterior derivative. In this case, ``\langle \cdot, \cdot \rangle`` denotes a chain-cochain pairing, which can be interpreted as integrating a cochain over a chain.

```math
\partial [v_0, v_1, \ldots, v_k] = \sum_{i=0}^{k} (-1)^i [v_0, \ldots, \hat{v}_i, \ldots, v_k]
```

The formula above displays how the boundary operator is computed in the discrete setting. The operator inputs a chain with a particular ordering of vertices, where edges may look like ``e = [v_0, v_1]`` and faces may look like ``\sigma = [v_0, v_1, v_2]``. Suppose the operator inputs the triangular face ``\sigma_0 = [v_0, v_1, v_2]`` surrounded by the edges ``e_0 = [v_0, v_1]``, ``e_1 = [v_1, v_2]``, and ``e_2 = [v_0, v_2]``. The boundary operator omits the ``i``-th vertex of the face at each step. The sign of each term is determined by whether ``i`` is even or odd. For this face ``\sigma_0``, the result is ``\partial [v_0, v_1, v_2] = [v_1, v_2] - [v_0, v_2] + [v_0, v_1] = e_1 - e_2 + e_0``. A sketch of the face with consistent edge orientations reproduces the same result.

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

#### Decapodes.jl

    :d₀ => dec_differential(0, sd) |> matmul
    :d₁ => dec_differential(1, sd) |> matmul
    :dual_d₀ || :d̃₀ => dec_dual_derivative(0, sd) |> matmul
    :dual_d₁ || :d̃₁ => dec_dual_derivative(1, sd) |> matmul

#### CombinatorialSpaces.jl

    d(s::HasDeltaSet, x::DualForm{n}) where n
    dual_derivative(n::Int, s::HasDeltaSet, args...) dual_derivative(::Type{Val{n}}, s::HasDeltaSet, args...) where n
    dual_derivative_nz(::Type{Val{0}}, s::AbstractDeltaDualComplex1D, x::Int) dual_derivative_nz(::Type{Val{0}}, s::AbstractDeltaDualComplex2D, x::Int) dual_derivative_nz(::Type{Val{1}}, s::AbstractDeltaDualComplex2D, x::Int)
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

The exterior derivative operator extends the idea of differentiation into differential geometry. It generalizes vector calculus concepts such as gradient, divergence, and curl to differential forms of various degrees. It is a central operator in many partial differential equations represented in the exterior calculus framework.

This operator linearly maps ``k``-forms to ``(k+1)``-forms (or ``k``-cochains to ``(k+1)``-cochains in the discrete setting), representing the form's rate of change on the mesh. In 3D space, the exterior derivative operator represents the gradient on a 0-form, the curl on a 1-form, and divergence on a 2-form. Applied to a 3-form in 3D space, this operator is always zero.

In the discrete setting, this operator matrix is determined by Generalized Stokes' Theorem, which states that integrating the exterior derivative of a form over a manifold is equivalent to integrating the form over the boundary of that manifold. Accordingly, the exterior derivative for a cochain is the transpose of the boundary matrix for the corresponding chain complex. Like the boundary operator, it is a sparse matrix with entries 0, +1, and -1. The computation of this operator is local, topological, and coordinate-free.

This operator determines the dimensional rate of change of a form based on the values of nearby forms. For example, imagine three cochains attached to 1-simplices that form a triangle. The exterior derivative of the resulting 1-cochain is a 2-cochain residing on the triangular face, computed by adding or subtracting the cochain values according to face and edge orientations.

Furthermore, the exterior derivative is nilpotent, meaning that applying it twice to the same form always gives zero, or ``d^2 = 0``. This property preserves geometric structure and mathematical conservation laws. For example, it encodes identities such as the fact that the divergence of a curl is zero and that the curl of a gradient is zero.

On a dual mesh, the discrete exterior derivative may change sign to account for the mesh's differing orientation. This correction appears in the Important Properties section.

Finally, this operator is natural with respect to discrete pullbacks, like the boundary.

### Important Properties

```math
d^k: \Omega_d^k(K) \to \Omega_d^{k+1}(K)
```

The discrete exterior derivative linearly maps discrete ``k``-forms on mesh ``K`` to discrete ``(k+1)``-forms on mesh ``K``.

```math
d(d\omega) = 0
```

The exterior derivative applied twice to the same form is always zero. This property retains the PDE's mathematical structure (geometric invariants, ``\nabla \times \nabla f = 0``, and ``\nabla \cdot \nabla \times f = 0``) without requiring any additional features.

```math
d(\omega \wedge \beta) = d\omega \wedge \beta + (-1)^{\deg(\omega)} \omega \wedge d\beta
```

The formula above displays the Leibniz Product Rule for the exterior derivative. The operator has a product-rule-like relationship with the wedge product.

```math
d = \partial^T
```

In the discrete setting, the exterior derivative matrix is computed as the transpose of the boundary operator. It inputs ``k``-forms and outputs the exterior derivative values of ``(k+1)``-forms.

```math
\widetilde{d}_{n-k} = (-1)^k (d_{k-1})^T
```

Due to the orientation changes in dual meshes, the exterior derivative applied to dual forms must be altered to account for differing orientation and orthogonality. The formula above displays this correction, where ``\widetilde{d}`` represents the dual exterior derivative operator.

```math
\int_{\Omega} d\omega = \int_{\partial \Omega} \omega
```

The formula above is Generalized Stokes' Theorem, which states that integrating the exterior derivative of a form over a manifold is equivalent to integrating the form over the boundary of that manifold. Thus, the exterior derivative generalizes the gradient, divergence, and curl, while Generalized Stokes' Theorem unifies Stokes' Theorem, Green's Theorem, and the Divergence Theorem into one expression.

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

#### Decapodes.jl

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

#### CombinatorialSpaces.jl

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

The wedge product operator is used to combine a ``k``-form and an ``l``-form into a higher-degree ``(k+l)``-form. As an example, imagine two vector-like 1-cochains. Their wedge product can be interpreted as the signed area formed by the two inputs, analogous to a determinant.

This operator is commutative when either ``k`` or ``l`` are even, meaning ``\beta \wedge \omega = \omega \wedge \beta``. However, if both forms have odd degrees, then the operator is anti-commutative, meaning ``\beta \wedge \omega = -\omega \wedge \beta``. It is also associative, where ``(\omega \wedge \beta) \wedge \gamma = \omega \wedge (\beta \wedge \gamma)``. In some DEC discretizations, including the one used in Decapodes, associativity is limited to closed forms. A closed form is a form whose exterior derivative is zero. Later operators that use this tool as a component also inherit this limitation.

Another limitation is that this discrete DEC wedge construction depends on the metric. In particular, its value is computed using support-volume measurements even though the smooth wedge product itself is topological.

Finally, the Leibniz Product Rule highlights an important relationship between the wedge product and exterior derivative. This rule is approximated using the operator formulas shown below.

### Important Properties

```math
\wedge : \Omega_d^k(K) \times \Omega_d^l(K) \to \Omega_d^{k+l}(K)
```

The discrete primal-primal wedge product inputs a discrete primal ``k``-form and discrete primal ``l``-form on mesh ``K`` and outputs a discrete ``(k+l)``-form on mesh ``K``. The discrete dual-dual wedge product and discrete primal-dual wedge product would be defined differently.

```math
\langle \omega^k \wedge \beta^l, \sigma^{k+l} \rangle :=
\frac{1}{(k+l)!}
\sum_{r \in S_{k+l+1}} \operatorname{sign}(\tau) \frac{|\sigma^{k+l} \cap \star v_{\tau(k)}|}{|\sigma^{k+l}|}
\langle \omega, [v_{\tau(0)}, \ldots, v_{\tau(k)}] \rangle
\langle \beta, [v_{\tau(k)}, \ldots, v_{\tau(k+l)}] \rangle
```

This is the metric formulation for the discrete primal-primal wedge product operator, where ``\omega`` is a discrete primal ``k``-form and ``\beta`` is a discrete primal ``l``-form. The left-hand side states that the evaluation of the wedge product of two forms on the appropriate face is equal to the right-hand side. The right-hand side states that for each permutation of the simplex ``\sigma^{k+l}`` (where there are ``(k+l+1)!`` orderings of the vertices), it finds the sign of that permutation, ``\operatorname{sign}(\tau)``, and multiplies it by a geometric weight, ``(|\sigma^{k+l} \cap \star v_{\tau(k)}| / |\sigma^{k+l}|)``. The weight determines the volume of the new simplex that is shared by the simplices touching the ``k``-th vertex of the current permutation, divided by the volume of the new simplex. Finally, it multiplies this weight by the evaluation of ``\omega`` on a ``k``-simplex formed from the initial vertices, ``\langle \omega, [v_{\tau(0)}, \ldots, v_{\tau(k)}] \rangle``, and the evaluation of ``\beta`` on an ``l``-simplex formed from the remaining vertices, ``\langle \beta, [v_{\tau(k)}, \ldots, v_{\tau(k+l)}] \rangle``. The contributions from each permutation are summed together and normalized by the factor ``(1/(k+l)!)``. As this formula relies on simplex magnitudes, it is metric.

```math
\langle \omega^k \wedge \beta^l, \sigma^{k+l} \rangle :=
\frac{1}{(k+l+1)!}
\sum_{r \in S_{k+l+1}} \operatorname{sign}(\tau)
\langle \omega, [v_{\tau(0)}, \ldots, v_{\tau(k)}] \rangle
\langle \beta, [v_{\tau(k)}, \ldots, v_{\tau(k+l)}] \rangle
```

The formula above is the topological primal-primal wedge product operator. It follows the same trend as the metric case. However, this one does not use a geometric weighting factor, and instead accounts for weight with the normalization ``(1/(k+l+1)!)``.

```math
\omega^k \wedge \beta^l = (-1)^{kl}(\beta^l \wedge \omega^k)
```

This formula expresses the graded-commutativity of the wedge product. If the form degrees ``k`` and ``l`` are both odd, switching the order of the forms swaps the sign of the resulting form. This reflects the fact that forms behave like antisymmetric tensors, so swapping the inputs that define a higher-degree form changes the orientation of the resulting area, volume, or higher-dimensional element. However, if either ``k`` or ``l`` is even, then swapping the terms does not change the sign.

```math
(\omega \wedge \beta) \wedge \gamma = \omega \wedge (\beta \wedge \gamma)
```

This formula displays the associative property of the wedge product operator in the smooth case. In the discrete case, this is only true for closed forms. The exterior derivative of a closed form is zero.

```math
d(\omega \wedge \beta) = d\omega \wedge \beta + (-1)^{\deg(\omega)} \omega \wedge d\beta
```

As displayed in the exterior derivative section, the Leibniz Product Rule defines the interaction between the exterior derivative and the wedge product operators.

```math
\omega \wedge \omega = 0
```

This formula shows the alternating property of the wedge product in the smooth case. Note that this vanishing property holds for odd-degree forms, by graded-commutativity.

### Citations

Discrete Exterior Calculus - Hirani\
Section 7.1-7.3 (pp. 71-76)

Discrete Exterior Calculus - Desbrun et. al.\
Section 8 (pp. 17-21)

## Hodge Star

### Signature

#### Decapodes.jl

    :⋆₀ => dec_hodge_star(0, sd, hodge=hodge) |> matmul
    :⋆₁ => dec_hodge_star(1, sd, hodge=hodge) |> matmul
    :⋆₂ => dec_hodge_star(2, sd, hodge=hodge) |> matmul
    :⋆₀⁻¹ => dec_inv_hodge_star(0, sd, hodge) |> matmul
    :⋆₁⁻¹ => dec_pair_inv_hodge(Val{1}, sd, hodge)
    :⋆₂⁻¹ => dec_inv_hodge_star(1, sd, hodge) |> matmul

#### CombinatorialSpaces.jl

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

The Hodge star operator maps ``k``-forms to their respective "orthogonal" ``(n-k)``-forms based on the [Riemannian metric](https://en.wikipedia.org/wiki/Riemannian_metric) (or inner product). Here, ``n`` is the number of dimensions in the space. It relates the primal mesh to the dual mesh by relating a chain to a corresponding dual chain. In 2D, primal vertices correspond to dual 2-cells, primal edges correspond to transverse dual edges, and primal faces correspond to dual vertices. It is also an important component in later operators.

The discrete Hodge star is computed using a formula that states that the original form evaluated on the original simplex, divided by the simplex's magnitude, is equivalent to the starred form evaluated on the orthogonal simplex, divided by the orthogonal simplex's magnitude. This formula ensures that the dual and primal meshes contain the same information. The operator is also an isomorphism, meaning that the Hodge star applied twice to a form returns the original form, up to a potential sign change depending on dimension.

The operator matrix can be computed in multiple ways but namely as the geometric Hodge or the diagonal Hodge. The geometric Hodge is calculated as a sparse, symmetric matrix with entries determined by the ratio between simplex and orthogonal simplex magnitudes. It pairs well with barycentric dual mesh centering. The matrix dimensions are ``|C_k| \times |C_k|``. The diagonal Hodge is a diagonal matrix with the same dimensions. It is faster to compute but less accurate. It pairs well with circumcentric dual mesh centering.

This operator is metric-dependent. It depends on a metric to determine angles and orthogonality. Thus, the result of this operator may differ even if simplex connectivity, cochain values, and orientations remain the same. For example, depending on whether the dual mesh is barycentric-centered or circumcentric-centered, the Hodge star may return different results. Furthermore, stretching the mesh changes simplex magnitudes on the primal and dual meshes, which in turn affects the output values of the Hodge star operator.

### Important Properties

```math
\star^k : \Omega_d^k(K) \to \Omega_d^{n-k}(\star K)
```

The discrete Hodge star linearly maps discrete ``k``-forms on mesh ``K`` to discrete ``(n-k)``-forms on dual mesh ``\star K``, where ``n`` is the number of dimensions of the space. The ``\star`` symbol acts on geometry and outputs its dual version.

```math
\langle \omega^k, \beta^k \rangle \, \mathrm{vol} = \omega^k \wedge \star \beta^k
```

The formula above displays the defining identity of the Hodge star in the smooth case. The inner product of two ``k``-dimensional forms, multiplied by the volume form, is equal to the wedge product between the first form and the Hodge star of the second form. The volume form encodes the metric and orientation, ensuring that the left-hand side has the same dimensions and sign as the right-hand side. Due to this definition, the Hodge star is inherently metric-dependent.

```math
\frac{1}{|\star \sigma^k|}
\langle \star \omega, \star \sigma^k \rangle \coloneqq
\frac{1}{|\sigma^k|} \langle \omega, \sigma^k \rangle
```

The formula above is the definition of the diagonal Hodge star for ``1 \le k \le n-1`` forms. It states that the Hodge star of a form integrated on the dual of a simplex and divided by the length of that dual simplex must be equal to the form evaluated on the primal simplex, divided by the length of that primal simplex. The dual form maintains the same physical information while transitioning to orthogonality. Forms that do not conform to the inequality have an additional term that changes sign to account for simplex orientation.

```math
\star\star \omega = (-1)^{k(n-k)} \omega
```

The Hodge star operator is an isomorphism. Using the Hodge star on a form twice returns it to its original state, with a change in sign depending on dual orientation.

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

#### Decapodes.jl

    :♭ᵈᵖ => dec_♭(sd)
    :♭♯ => ♭♯_mat(sd) |> matmul

#### CombinatorialSpaces.jl

    ♭(s::HasDeltaSet, X::DualVectorField)
    ♭(s::AbstractDeltaDualComplex2D, X::AbstractVector, ::DPPFlat)
    ♭(s::AbstractDeltaDualComplex2D, X::AbstractVector, ::PPFlat)
    ♭_mat(s::AbstractDeltaDualComplex2D, f::DPPFlat)
    ♭_mat(s::AbstractDeltaDualComplex2D, p2s, ::DPPFlat)
    ♭_mat(s::AbstractDeltaDualComplex2D, ::PPFlat)

### Description

The flat operator is a [musical isomorphism](https://en.wikipedia.org/wiki/Musical_isomorphism) that linearly maps vector fields to 1-forms. In the smooth case, it is the inverse of the sharp operator. In the discrete case, it has no exact inverse. It is often used to represent vector quantities, such as fluid flow, as forms.

The smooth definition of the flat operator depends on an inner product, making the operator metric-dependent. Because there are primal meshes, dual meshes, and multiple interpolation methods, there are eight different discrete flat operators. One example is the ``\flat_{dpp}`` operator, or dual-primal-primal operator. This operator inputs discrete dual vector fields, uses dual-primal interpolation, and outputs cochains on the primal mesh. For non-flat meshes, operators that input dual vector fields are often more appropriate. This operator is represented as a sparse matrix.

### Important Properties

```math
\flat : \mathfrak{X}_d(K) \to \Omega_d^1(K)
```

The discrete flat operator maps discrete vector fields, usually tangent vectors at primal vertices or circumcenters/barycenters, on mesh ``K`` to discrete 1-forms on ``K``.

```math
X^{\flat}(\cdot) = \langle X, \cdot \rangle
```

The formula defines the flat operator in the smooth case. Essentially, since ``X`` is given as a vector, the most natural way to convert it into a form is to plug it into an inner product. Note that this aligns with the use of the metric tensor to convert between vectors and forms.

```math
\langle X^{\flat_{dpp}}, \sigma^1 \rangle =
\sum_{\sigma^n \succ \sigma^1} \frac{|\star \sigma^1 \cap \sigma^n|}{|\star \sigma^1|}
X(\sigma^n) \cdot \tilde{\sigma^1}
```

The formula above displays the definition for the DPP flat operator, or the operator that inputs a vector field on the dual mesh, uses dual-primal interpolation, and outputs a cochain on the primal mesh. The LHS states that the flattened vector field evaluated on an edge is equivalent to the RHS. The RHS states that for every simplex that contains the evaluated edge ``(\sum_{\sigma^n \succ \sigma^1})``, determine how much of the evaluated edge's dual cell is within the simplex, find that magnitude, and then divide it by the magnitude of the edge's dual cell ``(|\star \sigma^1 \cap \sigma^n| / |\star \sigma^1|)``. Then, multiply this value by the dot product of the simplex's average vector with the unit vector of the evaluated edge ``(X(\sigma^n) \cdot \tilde{\sigma^1})``. Finally, sum all contributions.

### Citations

Discrete Exterior Calculus - Hirani\
Section 5.3-5.6 (pp. 46-54)

Discrete Exterior Calculus - Desbrun et. al.\
Section 7 (p. 16)

Notes on Discrete Exterior Calculus - Gillette\
Section 2.12 (pp. 14-15)

## Sharp

### Signature

#### Decapodes.jl

    :♯ᵖᵈ => dec_♯_pd(sd)
    :♯ᵖᵖ => dec_♯_pp(sd)
    :♯ᵈᵈ => dec_♯_dd(sd)

#### CombinatorialSpaces.jl

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

Because the smooth definition of the sharp depends on a metric, the sharp is a metric-dependent operator. There are four different discrete sharp operators to account for form inputs and vector outputs on the primal and dual meshes. It is represented as a sparse matrix.

### Important Properties

```math
\sharp : \Omega_d^1(K) \to \mathfrak{X}_d(K)
```

The discrete sharp operator maps discrete 1-forms on mesh ``K`` to discrete vector fields on mesh ``K``. The vectors are usually tangent to primal vertices or mesh circumcenters/barycenters.

```math
\langle \omega^{\sharp}, v \rangle = \omega(v)
```

The formula defines the sharp operator in the smooth case. The vector field ``\omega^{\sharp}`` is the unique vector field satisfying ``\langle \omega^{\sharp}, v \rangle = \omega(v)`` for every vector field ``v``. This relation holds at every point of a Riemannian manifold ``M``.

### Citations

Discrete Exterior Calculus - Hirani\
Chapter 5.7-5.8 (pp. 54-56)

Discrete Exterior Calculus - Desbrun et. al.\
Section 7 (p. 16)

## Codifferential

### Signature

#### Decapodes.jl

    :δ₁ => add_Codiff!(d, src, tgt)
    :δ₂ => add_Codiff!(d, src, tgt)

#### CombinatorialSpaces.jl

    δ(s::HasDeltaSet, x::SimplexForm{n}; kw...) where n
    δ(n::Int, s::HasDeltaSet, args...; kw...)
    δ(::Type{Val{n}}, s::HasDeltaSet; hodge::DiscreteHodge=GeometricHodge(), matrix_type::Type=SparseMatrixCSC{Float64}) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, form::AbstractVector; hodge::DiscreteHodge=GeometricHodge()) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, ::DiagonalHodge, args...) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, ::GeometricHodge, matrix_type) where n
    δ(::Type{Val{n}}, s::HasDeltaSet, ::GeometricHodge, form::AbstractVector) where n

### Description

The codifferential operator acts as the adjoint of the exterior derivative. It linearly maps ``k``-forms to ``(k-1)``-forms. This operator can represent divergence, and it is used in the algebraic definition of the Laplacian operator.

This operator is metric-dependent because it is defined using the Hodge star. This matches the smooth case, where it is the ``L^2`` adjoint of the exterior derivative and therefore requires a metric. It is represented as a sparse matrix created by combining Hodge star and exterior derivative matrices. Due to its simple definition, it can be computed by applying three matrix operators and accounting for orientation.

Finally, because the codifferential is defined using the exterior derivative, the smooth codifferential is nilpotent: applying it twice gives zero. However, due to inaccuracies in the discrete Hodge star, the discrete codifferential applied twice is only approximately zero.

### Important Properties

```math
\delta^k : \Omega_d^{k+1}(K) \to \Omega_d^k(K)
```

The discrete codifferential operator is a linear map that inputs discrete ``(k+1)``-forms on mesh ``K`` and outputs discrete ``k``-forms on ``K``.

```math
\delta f = 0
```

The formula above shows that the codifferential of a 0-form ``f`` is equal to zero. The codifferential operator takes a ``k``-form and outputs a ``(k-1)``-form. There is no such thing as a ``(-1)``-form, so the operation outputs zero.

```math
\delta \omega = (-1)^{n(k-1)+1} \star d \star \omega
```

The formula above displays the discrete definition of the codifferential operator on forms. It can be defined purely through Hodge star and exterior derivative operators. The ``(-1)^{n(k-1)+1}`` represents the orientation of the output. It is important to realize that if the codifferential is applied to a primal form, then the dual exterior derivative formula should be used in the definition, as ``\star \omega`` would produce a dual form.

```math
\langle d\omega, \beta \rangle = \langle \omega, \delta \beta \rangle
```

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

#### Decapodes.jl

    :i₁ => add_Inter_Prod_1D! or add_Inter_Prod_2D!(Val{1}, ...)
    :i₂ => add_Inter_Prod_2D!(Val{2}, ...)

    :ι₁₁ => interior_product_dd(Tuple{1,1}, sd)
    :ι₁₂ => interior_product_dd(Tuple{1,2}, sd)

#### CombinatorialSpaces.jl

    interior_product(s::HasDeltaSet, X♭::EForm, α::DualForm{n}; kw...) where n
    interior_product_flat(n::Int, s::HasDeltaSet, args...; kw...)
    interior_product_flat(::Type{Val{n}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...) where n
    interior_product_dd(::Type{Tuple{1,1}}, s::SimplicialSets.HasDeltaSet)
    interior_product_dd(::Type{Tuple{1,2}}, s::SimplicialSets.HasDeltaSet)

### Description

The interior product contracts a vector field with a ``k``-form, outputting a ``(k-1)``-form. It is an operator that "combines" forms and vector fields together by plugging the field into one of the form's inputs. In exterior algebra, forms are measurement tools that input vectors and output scalar values. The interior product makes use of this definition to merge vector fields and forms together.

In the smooth setting, the interior product importantly does not require a metric to be defined. However, using the algebraic definition, this tool becomes metric due to its reliance on the flat and Hodge star operators.

Because of the variant of the discrete DEC wedge product used in its definition, the Leibniz Product Rule for this operator only applies for closed forms due to limits in associativity. This operator is used in the Lie derivative, which faces this limitation as well.

Hirani notes the existence of an extrusion definition for this operator which is not covered in this description. However, it should be known that algebraic contraction is dual to the idea of extrusion in a geometric setting.

### Important Properties

```math
i_X : \Omega_d^k(K) \to \Omega_d^{k-1}(K)
```

The discrete interior product operator ``i_X`` maps ``k``-forms on discrete mesh ``K`` to ``(k-1)``-forms on ``K``, where ``X`` is a vector field.

```math
i_X \omega(X_1, \ldots, X_{k-1}) = \omega(X, X_1, \ldots, X_{k-1})
```

The above definition shows the algebraic definition of the interior product. This operator contracts a vector field with a ``k``-form, inserting it into one of the slots of the form. The output is a ``(k-1)``-form.

```math
i_X \omega = (-1)^{k(n-k)} \star (\star \omega \wedge X^{\flat})
```

The above identity displays the definition of the interior product operator composed of Hodge star, wedge product, and flat operators. The discrete operator can be created in the discrete setting with discrete forms of these operators. The structure of the operator, namely ``(\star \omega \wedge X^{\flat})``, represents inserting the contribution of the vector field ``X`` into the form ``\omega``. The exterior Hodge star then returns the result to its proper dimension, ``(k-1)``. Finally, the factor ``(-1)^{k(n-k)}`` accounts for changes in orientation due to the Hodge star and wedge product operators. Due to the use of the discrete wedge product in this identity, the Leibniz Product Rule for contraction only applies to closed forms, where ``d\omega = 0``, since the discrete wedge product is only associative on closed forms.

```math
i_X(f) = 0
```

The interior product of a 0-form ``f`` with a vector field ``X`` is always zero.

```math
i_X(\omega^k \wedge \beta^l) = (i_X \omega) \wedge \beta + (-1)^k \omega \wedge (i_X \beta)
```

The interior product operator follows the Leibniz Product Rule. For non-closed forms, it is important to note that the discrete wedge product is not associative and thus errors in this rule may appear.

### Citations

Discrete Exterior Calculus - Hirani\
Section 8.2-8.3 (pp. 79-83)

Discrete Exterior Calculus - Desbrun et. al.\
Section 10 (pp. 24-26)

## Lie Derivative

### Signature

#### Decapodes.jl

    :L₀ => add_Lie_1D!(Val{0}, ...)
    :L₁ => add_Lie_1D!(Val{1}, ...)
    :L₂ => add_Lie_2D!(Val{2}, ...)
    :ℒ₁ => ℒ_dd(Tuple{1,1}, sd)

#### CombinatorialSpaces.jl

    ℒ(s::HasDeltaSet, X♭::EForm, α::DualForm{n}; kw...) where n
    lie_derivative_flat(n::Int, s::HasDeltaSet, args...; kw...)
    lie_derivative_flat(::Type{Val{0}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    lie_derivative_flat(::Type{Val{1}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    lie_derivative_flat(::Type{Val{2}}, s::HasDeltaSet, X♭::AbstractVector, α::AbstractVector; kw...)
    ℒ_dd(::Type{Tuple{1,1}}, s::SimplicialSets.HasDeltaSet)

### Description

The Lie derivative operator inputs a ``k``-form and outputs a ``k``-form representing the rate of change of the form as it moves in the direction of a vector field. It is analogous to the directional derivative in vector calculus. A common application of the Lie derivative is advection.

The Lie derivative's algebraic definition is displayed in Cartan's Magic Formula, which uses exterior derivative and interior product operators. This smooth formula is used to define the discrete primal-dual Lie derivative, where ``X`` is a discrete vector field on the primal mesh and ``\omega`` is a dual ``k``-form.

However, this method may not satisfy the Leibniz Product Rule for Lie derivatives, as it uses the discrete wedge product operator within the interior product operator. This wedge product is only associative on closed forms.

Hirani notes an alternative flow-out definition in Section 8.5 of his thesis, Discrete Exterior Calculus, which does not suffer from this lack of associativity. However, that definition is not covered here.

### Important Properties

```math
\mathcal{L}_X : \Omega_d^k(K) \to \Omega_d^k(K)
```

The discrete Lie derivative ``\mathcal{L}_X`` maps ``k``-forms on discrete mesh ``K`` to ``k``-forms on ``K``, where ``X`` is a vector field.

```math
\mathcal{L}_X \omega = i_X(d\omega) + d(i_X \omega)
```

This formula displays Cartan's Magic Formula, or an algebraic definition of the Lie derivative using the interior product and the exterior derivative operators. This formula helps define the discrete primal-dual Lie derivative operator on a simplicial complex. Here, ``X`` is a discrete primal vector field and ``\omega`` is a dual ``p``-form.

```math
\mathcal{L}_X(\omega \wedge \beta) = \mathcal{L}_X(\omega) \wedge \beta + \omega \wedge \mathcal{L}_X(\beta)
```

This is the Leibniz Product Rule of the Lie derivative, which displays the Lie derivative's relationship with wedge product in the smooth case. In the discrete case, if the algebraic definition of the interior product is used, this relationship might not always be correct. This is because the algebraic interior product is defined using the wedge product. As the wedge product is only associative on closed forms (in the discrete setting), there may be errors in this property.

### Citations

Discrete Exterior Calculus - Hirani\
Section 8.4-8.5 (pp. 83-85)

Discrete Exterior Calculus - Desbrun et. al.\
Section 10 (pp. 25-26)

## Laplacian Operator

### Signature

#### Decapodes.jl

    :Δ₀ => add_De_Rham_1D!(Val{0}, ...)
    :Δ₁ => add_De_Rham_1D!(Val{1}, ...) or add_De_Rham_2D!(Val{1}, ...)
    :Δ₂ => add_De_Rham_2D!(Val{2}, ...)
    :Δᵈ₀ => Δᵈ(Val{0}, sd)
    :Δᵈ₁ => Δᵈ(Val{1}, sd)
    :Δ₀⁻¹ => dec_inv_lap_solver(Val{0}, sd)

#### CombinatorialSpaces.jl

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
    Δᵈ(::Type{Val{0}}, s::SimplicialSets.HasDeltaSet)
    Δᵈ(::Type{Val{0}}, s::SimplicialSets.HasDeltaSet2D)
    Δᵈ(::Type{Val{0}}, s::SimplicialSets.HasDeltaSet3D)
    Δᵈ(::Type{Val{1}}, s::SimplicialSets.HasDeltaSet)
    dec_Δ⁻¹(::Type{Val{0}}, s::AbstractGeometricMapSeries; scheme::AbstractSubdivisionScheme = BinarySubdivision(), steps = 3, cycles = 5, alg = cg, μ = 2)

### Description

The Laplacian operator, or the Laplace-Beltrami operator on 0-forms, generalizes the traditional Laplacian to curved surfaces. It linearly maps ``k``-forms to ``k``-forms. The ordinary Laplacian determines how much a point in a function deviates from the average of the points surrounding it. With the positive-semidefinite sign convention used in this code, a positive Laplacian indicates that the point lies below the local average, while a negative Laplacian indicates that it lies above the local average.

This operator is defined with the exterior derivative and codifferential operators. For 0-forms (scalars), this operator simplifies from ``\Delta = \delta d + d \delta`` to ``\delta d``, as the codifferential of a 0-form is always zero.

It is represented as a sparse matrix. Due to its reliance on the codifferential, and thus the Hodge star, it is metric-dependent.

### Important Properties

```math
\Delta^k : \Omega_d^k(K) \to \Omega_d^k(K)
```

The discrete Laplace-deRham operator, the general form of the Laplace-Beltrami operator, maps ``k``-forms on discrete mesh ``K`` to ``k``-forms on the same mesh.

```math
\Delta = d\delta + \delta d
```

This is the general Laplace-deRham operator, which maps ``k``-forms to dimensionally equivalent ``k``-forms in the smooth case. The Laplace-Beltrami operator, ``\Delta f = \delta d f``, is a special case for 0-forms where ``d(\delta f) = 0``, since the codifferential of a 0-form is always zero.

```math
\langle \Delta f, \sigma^0 \rangle = \frac{1}{|\star \sigma_0|} \sum_{\sigma^1 = [\sigma^0, v]} \frac{|\star \sigma^1|}{|\sigma^1|} (f(v) - f(\sigma^0))
```

This formula evaluates a 0-form ``f`` at a primal vertex ``\sigma^0`` on a well-oriented triangular mesh. The mesh does not need to be flat. It states that ``\langle \Delta f, \sigma^0 \rangle`` is computed by summing the contributions of all edges incident to ``\sigma^0``, weighting each contribution by the ratio ``|\star \sigma^1| / |\sigma^1|``, multiplying by the difference ``f(v) - f(\sigma^0)``, and then normalizing by the area of the corresponding dual cell ``|\star \sigma^0|``. This formula is equivalent to a different geometric Laplace-Beltrami formula that uses cotangents and indices.

### Citations

Discrete Exterior Calculus - Hirani\
Section 6.4 (pp. 68-70)

Discrete Exterior Calculus - Desbrun et. al.\
Section 10 (pp. 26-27)

## References

[Discrete Exterior Calculus](https://www.cs.jhu.edu/~misha/Fall09/Hirani03.pdf) - Hirani

[Discrete Exterior Calculus](https://arxiv.org/abs/math/0508341) - Desbrun et. al.

[Notes on Discrete Exterior Calculus](https://math.arizona.edu/~agillette/research/decNotes.pdf) - Gillette

[Discrete Differential Forms for Computational Modeling](https://geometry.caltech.edu/pubs/DKT05.pdf) - Desbrun et. al.
