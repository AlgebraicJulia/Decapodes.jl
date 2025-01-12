# Defining Operators

Decapodes comes with a predefined list of operators for the DEC. However, new
operators may be defined by either composition or by entries to a `generate`
function.

The `generate` function returns a lookup table for operators relevant for the
model. In this case, we define a model

``` @example DEC
function generate(sd, symbol; hodge=GeometricHodge())
    op = @match symbol begin
        :# => x -> 
