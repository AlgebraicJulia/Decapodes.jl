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

### Description

### Important Properties

### Citations


## Hodge Star

### Signature

### Description

### Important Properties

### Citations



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