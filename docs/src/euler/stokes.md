# Let's get Stoked!

Import the necessary packages and define a decapode for the Stokes dynamics.
We're using a discrete exterior calculus formulation of the [Stokes equations](https://en.wikipedia.org/wiki/Stokes_flow) in 2D.
This flow is more familiarly given by the equation

```math
\mu\Delta u - \nabla p + f=0
```

subject to the incompressibility constraint ``\nabla \cdot u = 0``.

```@example stoke
  using Catlab
  using CombinatorialSpaces
  using CombinatorialSpaces.ExteriorCalculus
  using CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form
  using ComponentArrays
  using StaticArrays
  using DiagrammaticEquations
  using DiagrammaticEquations.Deca
  using Decapodes
  using LinearAlgebra
  using CairoMakie
  import CairoMakie: wireframe, mesh, Figure, Axis

  StokesDynamics = @decapode begin
    (P)::Form0 ## Pressure.
    (v)::Form1 ## Velocity.
    (φ)::Constant
    (μ)::Constant

    ∂ₜ(v) == v̇
    ∂ₜ(P) == Ṗ
    
    v̇ == μ * Δ(v)-d₀(P) + φ
    Ṗ == ⋆₀⁻¹(dual_d₁(⋆₁(v)))
  end 
  to_graphviz(StokesDynamics)
```

Next we construct a sequence of finer and finer meshes
on a unit square.

```@example stoke
s0 = triangulated_grid(1,1,1/4,1/4*sqrt(3)/2,Point3D)
fs = repeated_subdivisions(3,s0,triforce_subdivision_map);
wireframe(dom(fs[1]).delta_set)
```

Next, we construct the operator 
instantiating the Stokes equations as a linear system
on ``\Omega_0\oplus \Omega_1`` and evaluate it on
the vector corresponding to a constant laminar flow
with constant pressure to show how the time derivatives
of ``v`` and ``P`` vanish in this case.

```@example stoke
s = s0
len(s,e) = sqrt(sum((s[:point][s[:∂v0][e]]- s[:point][s[:∂v1][e]]) .^2))
diam(s) = minimum(len(s,e) for e in edges(s))
ε = diam(s0)
bvs(s) = findall(x -> abs(x[1]) < ε || abs(x[1]-1) < ε || x[2] == 0 || abs(x[2]-1)< ε*sqrt(3)/2, s[:point])
ivs(s) = filter(x -> !(x in bvs(s)), 1:nv(s))
bes(s) = begin bvs_s = bvs(s) ; 
  findall(x -> s[:∂v0][x] ∈ bvs_s ||  s[:∂v1][x] ∈ bvs_s , parts(s,:E)) end
ies(s) = begin b = bes(s);  [e for e in edges(s) if !(e in b)] end
sd = dualize(s,Circumcenter())
gensim(StokesDynamics)
f! = evalsim(StokesDynamics)(sd,nothing)
#g! = (du,u) -> begin f!(du,u,(φ=zeros(ne(s)),μ=1),0);  end
```
This was a long block

```@example stoke
using LinearOperators, LinearSolve
g!(du::ComponentArray,u::ComponentArray,_,_) = begin 
  #f! = evalsim(StokesDynamics)(sd,nothing)
  f!(du,u,(φ=zeros(ne(s)),μ=1),missing)
#  du.v[bes(sd)] .= u.v[bes(sd)]; du.P[bvs(sd)] .= u.P[bvs(sd)]
  end
mask(u::ComponentArray) = u 
#ComponentArray(v=u.v[ies(sd)],P=u.P[ivs(sd)])
function g(u,a,b) 
  u = ComponentArray(v=u[1:ne(sd)],P=u[ne(sd)+1:end])
  v = similar(u)
  g!(v,u,a,b)
  v
end

function g!(du,u,_,_)
  uc = ComponentArray(v=u[1:ne(sd)],P=u[ne(sd)+1:end])
  duc = similar(uc)
  g!(duc,uc,nothing,nothing)
  du .= collect(duc)
end

opr = FunctionOperator(g,1:nv(sd)+ne(sd),1:nv(sd)+ne(sd),islinear=true)
nopr = LinearOperator(Float64,nv(sd)+ne(sd),nv(sd)+ne(sd),false,false,g!)
ω = eval_constant_primal_form(sd,SVector{3}([1.,1.,1.]))
u = ComponentArray(v=ω,P=ones(nv(sd)))
#up = opr(u,(φ=zeros(ne(s)),μ=1),0)
#norm(up.v[ies(sd)])
#norm(up.P[ivs(sd)])
#norm(up)

using SparseArrays
function matrix(f) 
  n,m = size(f)
  A = spzeros(n,m)
  for i in 1:n
    e = spzeros(n)
    e[i] = 1
    A[:,i] = f*e
  end
  A
end
#I think that A is faster than opr by a factor proportional to the width
A = matrix(opr)
nA = matrix(nopr)
b = rand(length(u))
b′ = ComponentVector(v=rand(ne(sd)),P=rand(nv(sd)))
```

Confirming that this really happened for some good reason,
we consider an iteration on random input data:

```@example stoke
ω′ = rand(ne(s))
P′ = rand(nv(s))
u′ = mask(ComponentArray(v=ω′,P=P′))
norm(g(u′))
```

Now let's try to solve the linear problem ``Lx=b`` where ``L`` is the operator ``opr`` we defined above.

```@example stoke
#A pretty random thing in the range of `opr`
b = opr*ComponentVector(v=ω,P=ones(nv(sd)))
prob = LinearProblem(opr,b)
sol = solve(prob)
norm(b)
t = ComponentVector(v=sol.u[1:ne(sd)],P=sol.u[ne(sd)+1:end])

norm(opr*t-b)
```
So the above mysteriosuly doesn't work; especially mysteriously considering tha this does work:
```@example stoke

b′ = collect(b)
prob′ = LinearProblem(A,b′)
sol′ = solve(prob′)
norm(A*sol′.u-b′)
```

Let's try solving with `opr` using a specified algorithm.

```@example stoke
sol3 = solve(prob,KrylovJL_GMRES())
t′ = ComponentVector(v=sol3.u[1:ne(sd)],P=sol3.u[ne(sd)+1:end])
norm(opr*sol3.u-b′)/norm(b′)
```

No luck. Let's switch to work directly with `Krylov.jl`. 

```@example stoke
using Krylov
b′ = nopr* collect(ComponentVector(v=ω,P=ones(nv(sd))))
sol4 = gmres(nopr,b′)
norm(nopr*sol4[1]-b′)/norm(b′)
```


