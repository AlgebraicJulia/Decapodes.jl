# Simple Equations

```@setup INFO
include(joinpath(@__DIR__, "..", "docinfo.jl"))
info = DocInfo.Info()
```

This tutorial shows how to use Decapodes to represent simple equations. These aren't using any of the Discrete Exterior Calculus or CombinatorialSpaces features of Decapodes. They just are a reference for how to build equations with the `@decapodes` macro and see how they are stored as ACSets.

```@example Harmonic
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
```

The harmonic oscillator can be written in Decapodes in at least three different ways.


```@example Harmonic
oscillator = @decapode begin
  X::Form0
  V::Form0

  ∂ₜ(X) == V
  ∂ₜ(V) == -k(X)
end
```

The default representation is a tabular output as an ACSet. The tables are `Var` for storing variables (`X`) and their types (`Form0`). `TVar` for identifying a subset of variables that are the tangent variables of the dynamics (`Ẋ`). The unary operators are stored in `Op1` and binary operators stored in `Op2`. If a table is empty, it doesn't get printed.

Even though a diagrammatic equation is like a graph, there are no edge tables, because the arity (number of inputs) and coarity (number of outputs) is baked into the operator definitions.

You can also see the output as a directed graph. The input arrows point to the state variables of the system and the output variables point from the tangent variables. You can see that I have done the differential degree reduction from  `x'' = -kx` by introducing a velocity term `v`. Decapodes has some support for derivatives in the visualization layer, so it knows that `dX/dt` should be called `Ẋ` and that `dẊ/dt` should be called `Ẋ̇`.

```@example Harmonic
to_graphviz(oscillator)
```

In the previous example, we viewed negation and transformation by `k` as operators. Notice that `k` appears as an edge in the graph and not as a vertex. You can also use a 2 argument function like multiplication (`*`). With a constant value for `k::Constant`. In this case you will see `k` enter the diagram as a vertex and multiplication with `*` as a binary operator.

```@example Harmonic
oscillator = @decapode begin
  X::Form0
  V::Form0

  k::Constant

  ∂ₜ(X) == V
  ∂ₜ(V) == -k*(X)
end
```

This gives you a different graphical representation as well. Now we have the cartesian product objects which represent a tupling of two values.

```@example Harmonic
to_graphviz(oscillator)
```

You can also represent negation as a multiplication by a literal `-1`.

```@example Harmonic
oscillator = @decapode begin
  X::Form0
  V::Form0

  k::Constant

  ∂ₜ(X) == V
  ∂ₜ(V) == -1*k*(X)
end
```

Notice that the type bubble for the literal one is `ΩL`. This means that it is a literal. The literal is also used as the variable name.

```@example Harmonic
infer_types!(oscillator)
to_graphviz(oscillator)
```

We can allow the material properties to vary over time by changing `Constant` to `Parameter`. This is how we tell the simulator that it needs to call `k(t)` at each time step to get the updated value for `k` or if it can just reuse that constant `k` from the initial time step.

```@example Harmonic
oscillator = @decapode begin
  X::Form0
  V::Form0

  k::Parameter

  ∂ₜ(X) == V
  ∂ₜ(V) == -1*k*(X)
end
```

```@example Harmonic
infer_types!(oscillator)
to_graphviz(oscillator)
```

Often you will have a linear material where you are scaling by a constant, and a nonlinear version of that material where that scaling is replaced by a generic nonlinear function. This is why we allow Decapodes to represent both of these types of equations.

```@example INFO
DocInfo.get_report(info) # hide
```
