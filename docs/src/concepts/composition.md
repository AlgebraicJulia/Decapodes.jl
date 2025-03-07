# Composition

``` @example DEC
using Catlab 
using DiagrammaticEquations 
using Decapodes 
```

```@raw html
<details open>
  <summary><i>Dependencies</i></summary>
  <pre><code>
  using Catlab
  using DiagrammaticEquations
  using Decapodes
  </code></pre>
</details>
```


Decapodes composition is formally known as an "operad algebra". That means that we don't have to encode our composition in a single undirected wiring diagram (UWD) and then apply it. Rather, we can define several UWDs, compose those, and then apply those. Of course, since the output of oapply is another Decapode, we could perform an intermediate oapply, if that is convenient. In this tutorial we will learn how to use this capability to define more complex models, and how this modularity allows us to swap out models.

Besides it being convenient to break apart large UWDs into component UWDs, this hierarchical composition can enforce rules on our physical quantities.

### Example: Klausmeier

We will study how the [Klausmeier model](@ref "Klausmeier") uses the operad
algebra to describe a phenomenon of vegetation patterns can be governed by the
interaction of two different dynamics.

```@example DEC
# See Klausmeier Equation 2.a
Hydrodynamics = @decapode begin
  (n,w)::DualForm0
  dX::Form1
  (a,ν)::Constant

  ∂ₜ(w) == a - w - w * n^2 + ν * L(dX, w)
end

# See Klausmeier Equation 2.b
Phytodynamics = @decapode begin
  (n,w)::DualForm0
  m::Constant

  ∂ₜ(n) == w * n^2 - m*n + Δ(n)
end
nothing # hide
```

However let's for a moment pretend we don't have a very good understanding of
the hydrodynamics. We'll create an incorrect model:
```@example DEC
WrongHydrodynamics = @decapode begin
    (n,w)::DualForm0
    dX::Form1
    ∂ₜ(w) == L(dX, w)
end
nothing # hide
```

Now that we have our two component models, we can specify a means of composing them via a composition pattern. This defines placeholder models which accept variables as inputs. If the variables between any two models are the same, then we have asserted that they are shared. For example, in this `@relation` diagram, we assert that there are two models, `phyto` and `hydro`, which each share the `N` and `W` variables.

```@example DEC
# Specify Composition
compose_klausmeier = @relation () begin
  phyto(N, W)
  hydro(N, W)
end

draw_composition(compose_klausmeier)
```

We apply our composition pattern by plugging in component Decapodes, and specifying which internal quantities to share along edges. Decapodes are formalized via the field of Applied Category Theory. A practical consequence here is that we can view a Decapode as a sort of computation graph.

```@example DEC
# Apply Composition
klausmeier_cospan = oapply(compose_klausmeier,
                           [Open(Phytodynamics, [:n, :w]),
                            Open(WrongHydrodynamics, [:n, :w])])
Klausmeier = apex(klausmeier_cospan)
to_graphviz(Klausmeier)
```

Let's pretend our understanding of the hydrodynamics model has improved. Since
the composition pattern has not changed, we can just substitute the
`WrongHydrodynamics` model out for the correct one.

```@example DEC
# Apply Composition
klausmeier_cospan = oapply(compose_klausmeier,
                           [Open(Phytodynamics, [:n, :w]),
                            Open(Hydrodynamics, [:n, :w])])
Klausmeier = apex(klausmeier_cospan)
to_graphviz(Klausmeier)
```
