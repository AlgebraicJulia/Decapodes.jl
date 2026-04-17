# ASCII and Vector Calculus Operators

Decapodes is able to present equations which closely resemble their
representation in mathematics because of Julia's native Unicode support. Most
editors (such as Vim or Emacs) or IDEs (such as VSCode) support Julia plug-ins
with built-in support for typing Unicode.

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
| avg₀₁    | avg_01     | average values stored on endpoints of edges   |

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

