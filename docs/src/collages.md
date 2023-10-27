# Collages: Formalizing Boundary and Initial Conditions

Representing physical equations via formal diagrams (i.e. a Decapode) solves many challenges in automatic simulation generation. In practice, however, setting initial and boundary conditions is a major part of the workflow, and this process itself must be formalized. We do so by implementing the notion of a *collage* of Decapodes. This is directly inspired by the [Diagrammatic view of differential equations in physics](https://arxiv.org/abs/2204.01843). See §4 for relevant theoretical discussion.

## Layman's Definition: Semantics

We can quickly introduce a *collage* as simply taking two Decapodes, specifying a *morphism* between the two (which point from quantities in one to quantities in the other), and interpreting this entire collection of arrows as itself a single large Decapode diagram. In contrast to the composition patterns used to specify composition between different physical components, we will not merge any quantities here. Rather, we interpret these new arrows as providing information on boundary and initial conditions. This information is passed directly to our simulation generator to help handle loading initial conditions and boundary conditions.

## The Workflow

### Goals
We have a few goals for this collage-based workflow:
- Make it quick and easy to try out different initial conditions
- Make it quick and easy to try out different boundary conditions
- Make it quick and easy to swap out different meshes
- Separate the specification of the "location" of boundaries from the "values" those boundaries should take
- Separate the specification of distributions from the properties of the mesh
- Formalize these representations for the sake of future-proofness and consistency

And of course, if this workflow is not ergonomic for your particular circumstances, you can hand-craft your initial and boundary conditions, while still using the other simulation automation that the Decapodes framework offers.

### Higher-Order Functions
To enforce a consistent workflow, we employ the programming technique of using higher-order functions. In particular, we will create a collage, and call `simulation_helper` on it. The simulation helper returns 3 things: a boundary condition **loader**, an initial condition **loader**, and a simulation **generator**.

The boundary condition **loader** accepts a dictionary which maps the names of physical quantities to functions which evaluate them on an arbitrary mesh. This loader then returns a **generator** which takes in the mesh that you want to run your simulation on, and any other regular constants and parameters that you want passed to your simulation. When you execute the generator on a particular mesh, you get a named tuple of parameters that is suitable to be used with the `ODEProblem` interface.

Similarly, the initial condition **loader** accepts a dictionary which maps the names of physical quantities to functions which evaluate them on an arbitrary mesh. This loader then returns a **generator** which takes in the mesh that you want to run your simulation on. When you execute the generator on a particular mesh, you get a state vector that is suitable to be used with the `ODEProblem` interface.

The simulation **generator** is the output of `gensim`, with operators that apply your boundary conditions included for you. Executing the output of `gensim` on a mesh returns a simulation function that is suitable to be used with the `ODEProblem` interface.
