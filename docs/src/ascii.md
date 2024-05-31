# ASCII and Vector Calculus Operators

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "docinfo.jl"))
info = DocInfo.Info()
```

Some users may have trouble entering unicode characters like ⋆ or ∂ in their development environment. So, we offer the following ASCII equivalents. Further, some users may like to use vector calculus symbols instead of exterior calculus symbols where possible. We offer support for such symbols as well.

## ASCII Equivalents

| Unicode  | ASCII      | Meaning                                       |
| -------  | -----      | -------                                       |
| ∂ₜ       | dt         | derivative w.r.t. time                        |
| ⋆        | star       | Hodge star, generalizing transpose            |
| Δ        | lapl       | laplacian                                     |
| ∧        | wedge      | wedge product, generalizing the cross product |
| ⋆⁻¹      | star\_inv  | Hodge star, generalizing transpose            |
| ∘(⋆,d,⋆) | div        | divergence, a.k.a. ∇⋅                         |

## Vector Calculus Equivalents

| Vec      | DEC              | How to enter                |
| -------- | ---------------- | --------------------------  |
| grad     | d                | grad                        |
| div      | ∘(⋆,d,⋆)         | div                         |
| curl     | ∘(d,⋆)           | curl                        |
| ∇        | d                | \nabla                      |
| ∇ᵈ       | ∘(⋆,d,⋆)         | \nabla \<tab\> \\^d \<tab\> |
| ∇x       | ∘(d,⋆)           | \nabla \<tab\> x            |
| adv(X,Y) | ∘(⋆,d,⋆)(X∧Y)    | adv                         |

```@example INFO
DocInfo.get_report(info) # hide
```
