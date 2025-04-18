# Types in Decapodes

This page describes the different types present in a Decapode. Most of them are derived from the DEC and are presented with short physical interpretations of their stored values.

## Types and What They Mean

| Type      | Description | Symbol in Graphviz  |
| :-------   | :-----       | -------: |
| `Form0`     | Data on primal vertices. Typically defined as a scalar function over the points on a mesh. | Ω₀|
| `Form1`     | Data on primal edges. Measures the strength of a vector field along each edge. | Ω₁|
| `Form2`     | Data on primal triangles. This is generally associated with a measure of flux through the area defined by the triangle. | Ω₂|
| `DualForm0` | Data on dual vertices. Each dual vertex lies at the center of a single primal triangle. The chosen center is usually chosen as either the Circumcenter (for close to equiangular triangles) or the Barycenter (for more irregular and obtuse triangles). | Ω̃₀|
| `DualForm1` | Data on dual edges. These edges are typically perpendicular to primal edges and as such can measure the flux across the primal edge. | Ω̃₁|
| `DualForm2` | Data on dual areas. Note that these areas are not necessarily triangles. | Ω̃₂|
| `Literal`   | A number like `1` or `3.14`. The user does not need to set this. | ΩL|
| `Constant`  | A value that is set by the user and will not change for the duration of the simulation. Typically this value is given as a scalar representing some physical constant. | ΩC|
| `Parameter` | A function provided by the user that takes in a time `t` and outputs a value. Like `Constant`, this is typically a scalar that changes over time. | ΩP|
| `infer` | A placeholder when the type of a certain variable is not yet or will not be inferred. These can usually be ignored but good typing helps with debugging and may improve performance. | Ω• |
