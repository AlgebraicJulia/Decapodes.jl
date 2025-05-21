# <img src='docs/src/logo/decapodes_logo.svg' width='128' alt='Georges de Prhawnm'> Decapodes.jl

[![Development Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://algebraicjulia.github.io/Decapodes.jl/dev)
[![Code Coverage](https://codecov.io/gh/AlgebraicJulia/Decapodes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AlgebraicJulia/Decapodes.jl)
[![CI/CD](https://github.com/AlgebraicJulia/Decapodes.jl/actions/workflows/julia_ci.yml/badge.svg)](https://github.com/AlgebraicJulia/Decapodes.jl/actions/workflows/julia_ci.yml)

Decapodes.jl is a framework for developing, composing, and simulating physical systems.

Decapodes.jl is the synthesis of Applied Category Theory (ACT) techniques for formalizing and composing physics equations, and Discrete Exterior Calculus (DEC) techniques for formalizing differential operators.
[CombinatorialSpaces.jl](https://algebraicjulia.github.io/CombinatorialSpaces.jl/dev/) hosts tools for discretizing space and defining DEC operators on simplicial complexes, and [DiagrammaticEquations.jl](https://github.com/AlgebraicJulia/DiagrammaticEquations.jl) hosts tooling for representing the equations as formal ACT diagrams. This repository combines these two packages, compiling diagrams down to simulatable code.

By combining the power of ACT and the DEC, we seek to improve the scientific computing workflow. Decapodes simulations are [hierarchically composable](https://algebraicjulia.github.io/Decapodes.jl/dev/bsh/budyko_sellers_halfar/), generalize over [any type of manifold](https://algebraicjulia.github.io/Decapodes.jl/dev/ice_dynamics/ice_dynamics/), and are [performant and accurate](https://www.cise.ufl.edu/~luke.morris/cism.html) with a declarative domain specific language (DSL) that is [human-readable](https://algebraicjulia.github.io/Decapodes.jl/dev/klausmeier/klausmeier/#Model-Representation).

<img src='https://algebraicjulia.github.io/Decapodes.jl/dev/navier_stokes/vort.gif' width='512' alt='Solution of the Vorticity-Formulation of the Incompressible Navier-Stokes Equations'>

## Citation

The ideas behind this package are described in the paper:

> Morris, Baas, Arias, Gatlin, Patterson, Fairbanks. Decapodes: A diagrammatic tool for representing, composing, and computing spatialized partial differential equations. *Journal of Computational Science* 81 (2024).
> [arXiv:2401.17432](https://arxiv.org/abs/2401.17432).
> [DOI:10.1016/j.jocs.2024.102345](https://doi.org/10.1016/j.jocs.2024.102345).
