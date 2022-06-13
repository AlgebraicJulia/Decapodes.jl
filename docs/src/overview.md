# Introduction to DECAPODEs

Discrete Exterior Calculus Applied to Partial and Ordinary Differential
Equations (DECAPODE) is a diagrammatic language used to express systems of
ordinary and partial differential equations. The DECAPODE provides a visual
framework for understanding the coupling between variables within a PDE or ODE
system, and a combinatorial data structure for working with them. Below, we
provide a high-level overview of how DECAPODEs can be generated and interpreted.

## Your First DECAPODE

We begin with the most basic DECAPODE, one which only includes a single
variable. In the DECAPODE graphical paradigm, nodes represent variables and
arrows represent operators which relate variables to each other. Since the
DECAPODE applies this diagrammatic language specifically to the Discrete
Exterior Calculus (DEC), variables are typed by the dimension and orientation
of the information they contain. So a variable of type `Form0{X}` will be the
0-dimensional data points on the space `X`, or in a discrete context, the
values defined on points of a mesh `X`. Similarly, `Form1{X}` will be values
stored on edges of the mesh, and `Form2{X}` will be values stored on the
surfaces of the mesh. Below, we provide a very simple DECAPODE with just a
single variable `C`. In this example, we also provide the necessary imports,
and define a convenience function for visualization in later examples.
```@example DEC
using Decapodes, Decapodes.Diagrams
using Catlab.Present, Catlab.Graphics

Variable = @decapode Decapodes2D begin
  C::Form0{X}
end;

draw_equation(decapode) = to_graphviz(decapode, node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"4.0"))

draw_equation(Variable)
```

The resulting diagram contains a single node, showing the single variable in
this system. We can then add a second variable:

```@example DEC
TwoVariables = @decapode Decapodes2D begin
  C::Form0{X}
  dC::Form1{X}
end;

draw_equation(TwoVariables)
```

And then can add some relationship between them. In this case, we make an
equation which states that `dC` is the discrete derivative of `C`:

```@example DEC
Equation = @decapode Decapodes2D begin
  C::Form0{X}
  dC::Form1{X}

  dC == d₀{X}(C)
end;

draw_equation(Equation)
```

Here, the two nodes represent the two variables, and the arrow between them
shows how they are related by the discrete derivative.

##  A Little More Complicated

Now that we've seen how to construct a simple equation, it's time to move on to
some actual PDE systems! One classic PDE example is the diffusion equation.
This equation states that the change of concentration at each point is
proportional to the laplacian of the concentration. One issue that we run into
here, though, is that there isn't a "proportionality" operator in the default
DECAPODEs syntax `Decapodes2D`. Thus, in this next example, we will first
extend the `Decapodes2D` syntax and then define the DECAPODE for diffusion.

```@example DEC
@present DiffusionQuantities <: Decapodes2D begin
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
end;

Diffusion = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
  # Diffusion equation
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end;

draw_equation(Diffusion)
```

The resulting DECAPODE shows the relationships between the three variables with
the triangle diagram. Note that automatic layout of these labels can result in
confusion as to which edge each label corresponds, but any confusion can be
resolved by referring back to the original `@decapode` definition.

## Bring in the Dynamics


