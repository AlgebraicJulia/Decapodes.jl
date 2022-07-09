# Poissuille Flow for Fluid Mechanics

When modeling a fluid flowing in pipe, one can ignore the multidimensional structure of the pipe and approximate the system as a 1 dimensional flow along the pipe. The noslip boundary condition and the geometry of the pipe enter a 1D equation in the form of a resistance term.

```@example Poiseuille
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using CombinatorialSpaces.DiscreteExteriorCalculus: ∧
using Catlab.Theories
using Catlab.Present
using Catlab.Graphics
using Catlab.Syntax
using Catlab.CategoricalAlgebra
using LinearAlgebra

using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules

# Julia community libraries

using CairoMakie
using Decapodes.Debug
using OrdinaryDiffEq


""" Decapodes1D
A schema which includes any homomorphisms that may be added by the @decapode
macro.

TODO: This should be chipped away at as more of this tooling takes advantage
of the Catlab GAT system
"""

@present Decapodes1D(FreeExtCalc1D) begin
  X::Space
  proj₁_⁰⁰₀::Hom(Form0(X)⊗Form0(X),Form0(X))
  proj₂_⁰⁰₀::Hom(Form0(X)⊗Form0(X),Form0(X))
  proj₁_⁰⁰₀⁺::Hom(Form0(X)⊕Form0(X),Form0(X))
  proj₂_⁰⁰₀⁺::Hom(Form0(X)⊕Form0(X),Form0(X))
  proj₁_⁰¹₀::Hom(Form0(X)⊗Form1(X),Form0(X))
  proj₂_⁰¹₁::Hom(Form0(X)⊗Form1(X),Form1(X))
  proj₁_⁰¹₀⁺::Hom(Form0(X)⊕Form1(X),Form0(X))
  proj₂_⁰¹₁⁺::Hom(Form0(X)⊕Form1(X),Form1(X))
  proj₁_⁰⁰̃₀::Hom(Form0(X)⊗DualForm0(X),Form0(X))
  proj₂_⁰⁰̃₀̃::Hom(Form0(X)⊗DualForm0(X),DualForm0(X))
  proj₁_⁰⁰̃₀⁺::Hom(Form0(X)⊕DualForm0(X),Form0(X))
  proj₂_⁰⁰̃₀̃⁺::Hom(Form0(X)⊕DualForm0(X),DualForm0(X))
  proj₁_⁰¹̃₀::Hom(Form0(X)⊗DualForm1(X),Form0(X))
  proj₂_⁰¹̃₁̃::Hom(Form0(X)⊗DualForm1(X),DualForm1(X))
  proj₁_⁰¹̃₀⁺::Hom(Form0(X)⊕DualForm1(X),Form0(X))
  proj₂_⁰¹̃₁̃⁺::Hom(Form0(X)⊕DualForm1(X),DualForm1(X))
  proj₁_¹⁰₁::Hom(Form1(X)⊗Form0(X),Form1(X))
  proj₂_¹⁰₀::Hom(Form1(X)⊗Form0(X),Form0(X))
  proj₁_¹⁰₁⁺::Hom(Form1(X)⊕Form0(X),Form1(X))
  proj₂_¹⁰₀⁺::Hom(Form1(X)⊕Form0(X),Form0(X))
  proj₁_¹¹₁::Hom(Form1(X)⊗Form1(X),Form1(X))
  proj₂_¹¹₁::Hom(Form1(X)⊗Form1(X),Form1(X))
  proj₁_¹¹₁⁺::Hom(Form1(X)⊕Form1(X),Form1(X))
  proj₂_¹¹₁⁺::Hom(Form1(X)⊕Form1(X),Form1(X))
  proj₁_¹⁰̃₁::Hom(Form1(X)⊗DualForm0(X),Form1(X))
  proj₂_¹⁰̃₀̃::Hom(Form1(X)⊗DualForm0(X),DualForm0(X))
  proj₁_¹⁰̃₁⁺::Hom(Form1(X)⊕DualForm0(X),Form1(X))
  proj₂_¹⁰̃₀̃⁺::Hom(Form1(X)⊕DualForm0(X),DualForm0(X))
  proj₁_¹¹̃₁::Hom(Form1(X)⊗DualForm1(X),Form1(X))
  proj₂_¹¹̃₁̃::Hom(Form1(X)⊗DualForm1(X),DualForm1(X))
  proj₁_¹¹̃₁⁺::Hom(Form1(X)⊕DualForm1(X),Form1(X))
  proj₂_¹¹̃₁̃⁺::Hom(Form1(X)⊕DualForm1(X),DualForm1(X))
  proj₁_⁰̃⁰₀̃::Hom(DualForm0(X)⊗Form0(X),DualForm0(X))
  proj₂_⁰̃⁰₀::Hom(DualForm0(X)⊗Form0(X),Form0(X))
  proj₁_⁰̃⁰₀̃⁺::Hom(DualForm0(X)⊕Form0(X),DualForm0(X))
  proj₂_⁰̃⁰₀⁺::Hom(DualForm0(X)⊕Form0(X),Form0(X))
  proj₁_⁰̃¹₀̃::Hom(DualForm0(X)⊗Form1(X),DualForm0(X))
  proj₂_⁰̃¹₁::Hom(DualForm0(X)⊗Form1(X),Form1(X))
  proj₁_⁰̃¹₀̃⁺::Hom(DualForm0(X)⊕Form1(X),DualForm0(X))
  proj₂_⁰̃¹₁⁺::Hom(DualForm0(X)⊕Form1(X),Form1(X))
  proj₁_⁰̃⁰̃₀̃::Hom(DualForm0(X)⊗DualForm0(X),DualForm0(X))
  proj₂_⁰̃⁰̃₀̃::Hom(DualForm0(X)⊗DualForm0(X),DualForm0(X))
  proj₁_⁰̃⁰̃₀̃⁺::Hom(DualForm0(X)⊕DualForm0(X),DualForm0(X))
  proj₂_⁰̃⁰̃₀̃⁺::Hom(DualForm0(X)⊕DualForm0(X),DualForm0(X))
  proj₁_⁰̃¹̃₀̃::Hom(DualForm0(X)⊗DualForm1(X),DualForm0(X))
  proj₂_⁰̃¹̃₁̃::Hom(DualForm0(X)⊗DualForm1(X),DualForm1(X))
  proj₁_⁰̃¹̃₀̃⁺::Hom(DualForm0(X)⊕DualForm1(X),DualForm0(X))
  proj₂_⁰̃¹̃₁̃⁺::Hom(DualForm0(X)⊕DualForm1(X),DualForm1(X))
  proj₁_¹̃⁰₁̃::Hom(DualForm1(X)⊗Form0(X),DualForm1(X))
  proj₂_¹̃⁰₀::Hom(DualForm1(X)⊗Form0(X),Form0(X))
  proj₁_¹̃⁰₁̃⁺::Hom(DualForm1(X)⊕Form0(X),DualForm1(X))
  proj₂_¹̃⁰₀⁺::Hom(DualForm1(X)⊕Form0(X),Form0(X))
  proj₁_¹̃¹₁̃::Hom(DualForm1(X)⊗Form1(X),DualForm1(X))
  proj₂_¹̃¹₁::Hom(DualForm1(X)⊗Form1(X),Form1(X))
  proj₁_¹̃¹₁̃⁺::Hom(DualForm1(X)⊕Form1(X),DualForm1(X))
  proj₂_¹̃¹₁⁺::Hom(DualForm1(X)⊕Form1(X),Form1(X))
  proj₁_¹̃⁰̃₁̃::Hom(DualForm1(X)⊗DualForm0(X),DualForm1(X))
  proj₂_¹̃⁰̃₀̃::Hom(DualForm1(X)⊗DualForm0(X),DualForm0(X))
  proj₁_¹̃⁰̃₁̃⁺::Hom(DualForm1(X)⊕DualForm0(X),DualForm1(X))
  proj₂_¹̃⁰̃₀̃⁺::Hom(DualForm1(X)⊕DualForm0(X),DualForm0(X))
  proj₁_¹̃¹̃₁̃::Hom(DualForm1(X)⊗DualForm1(X),DualForm1(X))
  proj₂_¹̃¹̃₁̃::Hom(DualForm1(X)⊗DualForm1(X),DualForm1(X))
  proj₁_¹̃¹̃₁̃⁺::Hom(DualForm1(X)⊕DualForm1(X),DualForm1(X))
  proj₂_¹̃¹̃₁̃⁺::Hom(DualForm1(X)⊕DualForm1(X),DualForm1(X))
  sum₀::Hom(Form0(X)⊗Form0(X),Form0(X))
  sum₁::Hom(Form1(X)⊗Form1(X),Form1(X))
  sum₀̃::Hom(DualForm0(X)⊗DualForm0(X),DualForm0(X))
  sum₁̃::Hom(DualForm1(X)⊗DualForm1(X),DualForm1(X))
end
```

# Creating the Poiseuille Equations

The first step is to present an extension of the generic Decapodes1D presentation with specific named linear operators for the viscosity effect and the drag effect. This is purely syntactic, we will add the corresponding matrices later. 

The `@decapode` macro creates the data structure representing the equations of Poiseuille flow. The first block declares variables, the second block defines intermediate terms and the last block is the core equation.

```@example Poiseuille
@present Poiseuille <: Decapodes1D begin
  (R, μ̃)::Hom(Form1(X), Form1(X))
  # μ̃ = negative viscosity per unit area
  # R = drag of pipe boundary
end;

Poise = @decapode Poiseuille begin
  (∇P)::Form1{X}
  (q, q̇, Δq)::Form1{X}
  P::Form0{X}

  # Laplacian of q for the viscous effect
  Δq == d₀{X}(⋆₀⁻¹{X}(dual_d₀{X}(⋆₁{X}(q))))
  # Gradient of P for the pressure driving force
  ∇P == d₀{X}(P)
  # definition of time derivative of q
  ∂ₜ{Form1{X}}(q) == q̇

  # the core equation
  q̇ == sum₁(sum₁(μ̃(Δq), ∇P),R(q))
end;
```

# Defining the Semantics

In order to solve our equations, we will need numerical linear operators that give meaning to our symbolic operators. The operator funcs code below assigns the necessary matrices as definitions for the symbols. In order to define the viscosity effect correctly we have to identify boundary edges and apply a mask. This is because the DEC has discrete dual cells at the boundaries that need to be handled specially for the viscosity term. We found empirically that if you allow nonzero viscosity at the boundary edges, the flows at the boundaries will be incorrect. 

```@example Poiseuille
"""    boundary_edges(ds)

Compute the edges of a 1D simplicial set that are either incident to in-degree 1 or out-degree 1 nodes.
For a graph, these are boundary vertices meaning leaf nodes. For our pipeflow problems,
these are the edges where material can enter the pipe network.
"""
function boundary_edges(ds)
  out_degree(x) = length(incident(ds, x, :∂v1))
  in_degree(x) = length(incident(ds, x, :∂v0))
  bpoints = findall(x -> out_degree(x) == 0 || in_degree(x) == 0, 1:nv(ds))
  sedges = vcat(incident(ds,bpoints,:∂v0)...)
  tedges = vcat(incident(ds,bpoints,:∂v1)...)
  bedges = collect(union(sedges, tedges))
  return bedges
end

"""    mask_boundary_edges(ds)

Provides the `boundary_edges(ds)` as a vector of 0/1 entries to use as a mask.
"""
function mask_boundary_edges(ds)
  D = ones(Int, ne(ds))
  D[boundary_edges(ds)] .= 0
  return D
end


opbind(f, T) = Dict(:operator=>f, :type=>T)

function create_funcs(ds, hodge=DiagonalHodge())
  funcs = Dict{Symbol, Dict}()
  funcs[:⋆₁] = opbind(⋆(Val{1}, ds, hodge=hodge), MatrixFunc())
  funcs[:⋆₁] = opbind(⋆(Val{1}, ds, hodge=hodge), MatrixFunc())
  funcs[:⋆₀] = opbind(⋆(Val{0}, ds, hodge=hodge), MatrixFunc())
  funcs[:⋆₀⁻¹] = opbind(inv(⋆(Val{0}, ds, hodge=hodge)), MatrixFunc())
  funcs[:⋆₁⁻¹] = opbind(inv(⋆(Val{1}, ds, hodge=hodge)), MatrixFunc())
  funcs[:d₀] = opbind(d(Val{0}, ds), MatrixFunc())
  funcs[:dual_d₀] = opbind(dual_derivative(Val{0}, ds), MatrixFunc());
  funcs[:sum₁] = opbind((x′, x, y)->(x′ .= x .+ y), InPlaceFunc())
  funcs[:∧₀₁] = opbind((r, c, v) -> r .= -∧(Tuple{0,1}, ds, c, v), InPlaceFunc())
  return funcs
end

##
function create_funcs(ds, operators, boundaries, hodge=DiagonalHodge())
  funcs = create_funcs(ds, hodge)
  merge!(funcs, operators, boundaries)
  return funcs
end

function operator_funcs(ds)
  F = Dict(
    :μ̃ => opbind(0.5 *  Diagonal(mask_boundary_edges(ds)), MatrixFunc()),
    :R => opbind(-0.1 * I(ne(ds)), MatrixFunc()),
    :¬ => opbind(-I(ne(ds)), MatrixFunc()),
    :k => opbind(1.0 * I(nv(ds)), MatrixFunc()))
  B = Dict(
    :∂ρ => opbind((ρᵇ, ρ) -> begin ρᵇ .= ρ; ρᵇ[1] = 0; ρᵇ[end] = 0; return ρᵇ end, InPlaceFunc())
  )
  create_funcs(ds, F, B)
end
form2dim = Dict(:Scalar => x->1,
                :Form0 => nv,
                :Form1 => ne,
                :DualForm1 => nv,
                :DualForm0 => ne)

```

## A Single Pipe Segment

We create a mesh with one pipe segment to see if we get the right answer. This simulation can be validated with the Poiseuille equation for a single pipe. First we create the mesh.

```@example Poiseuille
Point3D = Point3{Float64}
s = EmbeddedDeltaSet1D{Bool,Point3D}()
add_vertices!(s, 2, point=[Point3D(-1, 0, 0), Point3D(+1, 0, 0)])
add_edge!(s, 1, 2, edge_orientation=true)

ds = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
subdivide_duals!(ds, Circumcenter())
ds
```

Then we solve the equations.

```@example Poiseuille
funcs = operator_funcs(ds)
func, code = gen_sim(diag2dwd(Poise), funcs, ds; autodiff=false, form2dim=form2dim, params=[:P]);
prob = ODEProblem(func, [2.], (0.0, 10000.0), [1.,11.])
sol = solve(prob, Tsit5(); progress=true);
sol.u
```

## A Linear Pipe with Multiple Segments

We then move on to a linear sequence of pipe segments. You can visualize this as the discretization of a single long pipe into `n` segments. First we define the mesh:

```@example Poiseuille
function linear_pipe(n::Int)
  s = EmbeddedDeltaSet1D{Bool,Point3D}()
  add_vertices!(s, n, point=[Point3D(i, 0, 0) for i in 1:n])
  add_edges!(s, 1:n-1, 2:n, edge_orientation=true)
  orient!(s)
  ds = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
  subdivide_duals!(ds, Circumcenter())
  funcs = operator_funcs(ds)
  func, _ = gen_sim(diag2dwd(Poise), funcs, ds; autodiff=false, form2dim=form2dim, params=[:P])
  return ds, func, funcs
end

ds, func, funcs = linear_pipe(10)
ds
```

Then we solve the equation. Notice that the equilibrium flow is constant down the length of the pipe. This must be true because of conservation of mass. The segments are all the same length and the total flow in must equal the total flow out of each segment.

```@example Poiseuille
prob = ODEProblem(func, [5,3,4,2,5,2,8,4,3], (0.0, 10000.0), [10. *i for i in 1:10])
sol = solve(prob, Tsit5(); progress=true);
sol.u
```

## A Distribution Network

To model a distribution network such as residential drinking water or natural gas, we will build a binary tree of pipes that at each junction have a bifurcation into two pipes. We expect that the flow will be divided by two at each level of the tree. First we make the mesh.

```@example Poiseuille
function binary_pipe(depth::Int)
  s = EmbeddedDeltaSet1D{Bool,Point3D}()
  add_vertex!(s, point=Point3D(0, 0, 0))
  for n in 1:depth
    for prev_v in vertices(s)[end-2^(n-1)+1:end]
      x, y, _ = s[:point][prev_v]
      vs = add_vertices!(s, 2, point=[Point3D(sgn*3^0.5 + x, y+1, 0)
                                 for sgn in [1,-1]])
      add_edges!(s, vs, [prev_v,prev_v], edge_orientation=true)
    end
  end
  v = add_vertex!(s, point=Point3D(3^0.5, -1, 0))
  add_edge!(s, 1, v, edge_orientation=true)
  orient!(s)
  ds = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
  subdivide_duals!(ds, Circumcenter())
  funcs = operator_funcs(ds)
  func, _ = gen_sim(diag2dwd(Poise), funcs, ds; autodiff=false, form2dim=form2dim, params=[:P])
  return ds, func, funcs
end
ds, func, funcs = binary_pipe(2);
ds
```

Then we solve the equations.

```@example Poiseuille
prob = ODEProblem(func,
                 [5. for _ in 1:ne(ds)],
                 (0.0, 10000.0),
                 Float64[2^(7-p[2]) for p in ds[:point]])

sol = solve(prob, Tsit5(); progress=true);
sol.u
```

## Multiphysics

Decapodes really shines when you want to extend or refine your physics. We will change our physics by adding in a term for density of the material and the corresponding changes in pressure. This is not the only formulation for including a dynamic pressure effect into this system. If you can think of a better way to include this effect, we invite you to try it as an exercise!

Because the pressure is no longer being supplied as a parameter of the system controlled by the operators, we need to introduce a density term and a boundary condition for that density. In this system you can think of forcing a prescribed amount of material per unit time through the openings of the pipe and allowing the flow (q) and the pressure (P) to fluctuate. Before we were enforcing a fixed pressure gradient and and letting the flow fluctuate to achieve equilibrium. In the prior model, we were not accounting for the amount of material that had to flow in order to achieve that (flow, pressure) combination.


The Decapode can be visualized with graphviz, note that the boundary conditions are explicitly represented in the Decapode as operators that implement a masking operation. This is not consistent with the Diagrammatic Equations in Physics paper [PBHF22]. This approach is more directly tied to the computational method and will eventually be replaced with one based on morphisms of diagrams.

```@example Poiseuille
@present Poiseuille <: Decapodes1D begin
  (R, μ̃, ¬)::Hom(Form1(X), Form1(X))
  k::Hom(Form0(X), Form0(X))
  # μ̃ = negative viscosity per unit area
  # R = drag of pipe boundary
  # k = pressure as a function of density
  # boundary conditions
  ∂ρ::Hom(Form0(X), Form0(X))
end;

Poise = @decapode Poiseuille begin
  (∇P)::Form1{X}
  (q, q̇, Δq)::Form1{X}
  (P, ρ, ρ̇)::Form0{X}

  # Poiseuille Flow
  Δq == d₀{X}(⋆₀⁻¹{X}(dual_d₀{X}(⋆₁{X}(q))))
  ∂ₜ{Form1{X}}(q) == q̇
  ∇P == d₀{X}(P)
  q̇ == sum₁(sum₁(μ̃(Δq), ¬(∇P)),R(q))
  
  # Pressure/Density Coupling
  P == k(ρ)
  ∂ₜ{Form0{X}}(ρ) == ρ̇
  ρ̇ == ⋆₀⁻¹{X}(dual_d₀{X}(⋆₁{X}(∧₀₁{X}(ρ,q)))) # advection
  
  # Boundary conditions
  ρᵇ::Form0{X}
  ∂ρ(ρ̇) == ρᵇ
end;

to_graphviz(Poise, node_labels=true, prog="neato", node_attrs=Dict(:shape=>"oval"))
```

Then we can create the mesh and solve the equation.

```@example Poiseuille
function linear_pipe(n::Int)
  s = EmbeddedDeltaSet1D{Bool,Point3D}()
  add_vertices!(s, n, point=[Point3D(i, 0, 0) for i in 1:n])
  add_edges!(s, 1:n-1, 2:n, edge_orientation=true)
  orient!(s)
  ds = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
  subdivide_duals!(ds, Circumcenter())
  funcs = operator_funcs(ds)
  func, _ = gen_sim(diag2dwd(Poise, in_vars=[:q, :ρ]), funcs, ds; autodiff=false, form2dim=form2dim, params=[:P])
  return ds, func, funcs
end

ds, func, funcs = linear_pipe(10)

prob = ODEProblem(func, [5,3,4,2,5,2,3,4,3, 10,9,8,7,6,5,5,5,5,5], (0.0, 10000.0), [10. *i for i in 1:10])
sol = solve(prob, Tsit5(); progress=true);
sol.u
```

Notice that the solution contains both a vector of flows and a vector of pressures.
