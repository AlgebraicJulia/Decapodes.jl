# Types in Decapodes

This page describes the different kinds of types present in a Decapode. Most of them are derived from the DEC and are presented with short physical interpretations of their stored values.

## Types and What They Mean

| Type      | Description | Symbol in Graphviz  |
| :-------   | :-----       | -------: |
| `Form0`     | Information on primal vertices. Typically defined as a scalar function over the points on a mesh. | Ω₀|
| `Form1`     | Information on primal edges. Essentially measuring the strength of a vector field along the edge.| Ω₁|
| `Form2`     | Information on primal triangles. This is generally associated with a measure of flux through the area defined by the triangle. | Ω₂|
| `DualForm0` | Information on dual vertices. Each dual vertex lies at the center of a single primal triangle. The chosen center is usually chosen as either the Circumcenter (for close to equiangular triangles) or the Barycenter (for more irregular and obtuse triangles).| Ω̃₀|
| `DualForm1` | Information on dual edges. These edges are typically perpendicular to primal edges and as such can measure the flux across the primal edge.| Ω̃₁|
| `DualForm2` | Information on dual areas. Note that these areas are not necessarily triangles. | Ω̃₂|
| `Literal`   | A number like `1` or `pi`. This type is set automatically by DECAPODES.| ΩL|
| `Constant`  | A value that is set by the user and will not change for the duration of the simulation. Typically this value is given as a scalar representing some physical constant. | ΩC|
| `Parameter` | A function provided by the user that takes in a time `t` and outputs a value. Like `Constant` this is typically a scalar that changes over time.| ΩP|
| `infer` | A type specifying that DECAPODES can not figure out the type of a certain variable. These can usually be ignored but good typing helps with debugging and may improve performance. | Ω• |

## Using Primal vs Dual

In the DEC there is a division between information stores in primal and dual forms. 

<!-- TODO: Explain! --->