# Poissuille Flow for Fluid Mechanics

When modeling a fluid flowing in a pipe, one can ignore the multidimensional structure of the pipe and approximate the system as a 1 dimensional flow along the pipe. The noslip boundary condition and the geometry of the pipe enter a 1D equation in the form of a resistance term.

```@example Poiseuille
using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using CombinatorialSpaces.DiscreteExteriorCalculus: ∧
using Decapodes

# Julia community libraries
using CairoMakie
using GeometryBasics
using LinearAlgebra
using OrdinaryDiffEq
```

# Creating the Poiseuille Equations

The `@decapode` macro creates the data structure representing the equations of Poiseuille flow. The first block declares variables, the second block defines intermediate terms and the last block is the core equation.

```@example Poiseuille
# μ̃ = negative viscosity per unit area
# R = drag of pipe boundary

Poise = @decapode begin
  P::Form0
  q::Form1
  (R, μ̃ )::Constant

  # Laplacian of q for the viscous effect
  Δq == Δ(q)
  # Gradient of P for the pressure driving force
  ∇P == d(P)

  # Definition of the time derivative of q
  ∂ₜ(q) == q̇

  # The core equation
  q̇ == μ̃  * ∂q(Δq) + ∇P + R * q
end

Poise = expand_operators(Poise)
infer_types!(Poise, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(Poise, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(Poise)
```

# Defining the Semantics

In order to solve our equations, we will need numerical linear operators that give meaning to our symbolic operators. The `generate` function below assigns the necessary matrices as definitions for the symbols. In order to define the viscosity effect correctly we have to identify boundary edges and apply a mask. This is because the DEC has discrete dual cells at the boundaries that need to be handled specially for the viscosity term. We found empirically that if you allow nonzero viscosity at the boundary edges, the flows at the boundaries will be incorrect. 

```@example Poiseuille
using MLStyle
include("../../examples/boundary_helpers.jl")

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :∂q => x -> begin
      x[boundary_edges(sd)] .= 0
      x
    end
    :∧₀₁ => (x,y) -> begin
      ∧(Tuple{(0,1)}, sd, x,y)
    end
    :∂ρ => ρ -> begin
      ρ[1] = 0
      ρ[end] = 0
      ρ
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end
```

## A Single Pipe Segment

We create a mesh with one pipe segment to see if we get the right answer. This simulation can be validated with the Poiseuille equation for a single pipe. First we create the mesh.

```@example Poiseuille
Point3D = Point3{Float64}
s = EmbeddedDeltaSet1D{Bool,Point3D}()
add_vertices!(s, 2, point=[Point3D(-1, 0, 0), Point3D(+1, 0, 0)])
add_edge!(s, 1, 2, edge_orientation=true)

sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())
sd
```

Then we solve the equations.

```@example Poiseuille
using MultiScaleArrays
sim = eval(gensim(Poise))
fₘ = sim(sd, generate)
q = [2.0]
P = [10.0, 5.0]
u = construct(PhysicsState, [VectorForm(q), VectorForm(P)], Float64[], [:q, :P])
params = (k = -0.01, μ̃ = 0.5, R=0.005)
prob = ODEProblem(fₘ, u, (0.0, 10000.0), params)
sol = solve(prob, Tsit5())
sol.u
```

## A Linear Pipe with Multiple Segments

We then move on to a linear sequence of pipe segments. You can visualize this as the discretization of a single long pipe into `n` segments. First we define the mesh:

```@example Poiseuille
function linear_pipe(n::Int)
  s = EmbeddedDeltaSet1D{Bool,Point3D}()
  add_vertices!(s, n, point=[Point3D(i, 0, 0) for i in 1:n])
  add_edges!(s, 1:n-1, 2:n, edge_orientation=true)
  sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
  subdivide_duals!(sd, Circumcenter())
  sd
end

sd = linear_pipe(10)
true # hide
```

Then we solve the equation. Notice that the equilibrium flow is constant down the length of the pipe. This must be true because of conservation of mass. The segments are all the same length and the total flow in must equal the total flow out of each segment.

Note that we do not generate new simulation code for Poiseuille flow with `gensim` again. We provide our new mesh so that our discrete differential operators can be instantiated.

```@example Poiseuille
fₘ = sim(sd, generate)
P = [9,8,7,6,5,4,3,2,1,0]
q = [5,3,4,2,5,2,8,4,3]
u = construct(PhysicsState, [VectorForm(q), VectorForm(P)], Float64[], [:q, :P])
params = (k = -0.01, μ̃ = 0.5, R=0.005)
prob = ODEProblem(fₘ, u, (0.0, 10000.0), params)
sol = solve(prob, Tsit5());
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
  sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
  subdivide_duals!(sd, Circumcenter())
  sd
end
sd = binary_pipe(2)
true # hide
```

Then we solve the equations.

```@example Poiseuille
fₘ = sim(sd, generate)
P = collect(1.0:nv(sd))
q = fill(5.0, ne(sd))
u = construct(PhysicsState, [VectorForm(q), VectorForm(P)], Float64[], [:q, :P])
params = (k = -0.01, μ̃ = 0.5, R=0.005)
prob = ODEProblem(fₘ, u, (0.0, 10000.0), params)
sol = solve(prob, Tsit5())
sol.u
```

## Multiphysics

Decapodes really shines when you want to extend or refine your physics. We will change our physics by adding in a term for density of the material and the corresponding changes in pressure. This is not the only formulation for including a dynamic pressure effect into this system. If you can think of a better way to include this effect, we invite you to try it as an exercise!

Because the pressure is no longer being supplied as a parameter of the system controlled by the operators, we need to introduce a density term and a boundary condition for that density. In this system you can think of forcing a prescribed amount of material per unit time through the openings of the pipe and allowing the flow (q) and the pressure (P) to fluctuate. Before we were enforcing a fixed pressure gradient and and letting the flow fluctuate to achieve equilibrium. In the prior model, we were not accounting for the amount of material that had to flow in order to achieve that (flow, pressure) combination.


The Decapode can be visualized with graphviz, note that the boundary conditions are explicitly represented in the Decapode as operators that implement a masking operation. This is not consistent with the Diagrammatic Equations in Physics paper [PBHF22]. This approach is more directly tied to the computational method and will eventually be replaced with one based on morphisms of diagrams.

```@example Poiseuille
# μ̃ = negative viscosity per unit area
# R = drag of pipe boundary
# k = pressure as a function of density
Poise = @decapode begin
  q::Form1
  (P, ρ)::Form0
  (k, R, μ̃ )::Constant

  # Poiseuille Flow
  Δq == ∘(d, ⋆, d, ⋆)(q)
  ∂ₜ(q) == q̇
  ∇P == d(P)
  q̇ == μ̃ * ∂q(Δq) - ∇P + R * q
  
  # Pressure/Density Coupling
  P == k * ρ
  ∂ₜ(ρ) == ρ̇
  ρ_up == ∘(⋆, d, ⋆)(-1 * ∧₀₁(ρ,q)) # advection
  
  # Boundary conditions
  ρ̇ == ∂ρ(ρ_up)
end

Poise = expand_operators(Poise)
infer_types!(Poise, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(Poise, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(Poise)
```

Then we can create the mesh and solve the equation.

```@example Poiseuille
# Create mesh and subdivide it.
function linear_pipe(n::Int)
  s = EmbeddedDeltaSet1D{Bool,Point3D}()
  add_vertices!(s, n, point=[Point3D(i, 0, 0) for i in 1:n])
  add_edges!(s, 1:n-1, 2:n, edge_orientation=true)
  sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point3D}(s)
  subdivide_duals!(sd, Circumcenter())
end

sd = linear_pipe(20)

sim = gensim(Poise)
func = sim(sd, generate)

q = [5,3,4,2,5,2,3,4,3, 10,9,8,7,6,5,5,5,5,5]
ρ = [5,3,4,2,5,2,3,4,3, 10,9,8,7,6,5,5,5,5]
u = construct(PhysicsState, [VectorForm(q), VectorForm(ρ)], Float64[], [:q, :ρ])
params = (k = -0.01, μ̃ = 0.5, R=0.005)

prob = ODEProblem(func, u, (0.0, 10000.0), params)
sol = solve(prob, Tsit5())
sol.u
```

Notice that the solution contains both a vector of flows and a vector of pressures.
