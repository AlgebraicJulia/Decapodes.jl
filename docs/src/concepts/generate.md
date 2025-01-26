``` @example DEC
using Decapodes
```

## Custom Operators

Decapodes.jl already defines a suite of operators from the Discrete Exterior Calculus. However it is often the case that an implementation requires custom operators. Sometimes, this is just matter of building operators through composition. However Decapodes accepts a lookup table of functions which are included when parsing a Decapodes expression. This allows for new operators with their own symbols to be defined.

On this page, we will give an overview of the Decapodes `generate` function. 

### The `generate` function

The `gensim` function optionally accepts a callable object like a function to act as a lookup table for new operators. In general practice, this function is called `generate`, but this is not necessary. It just requires as arguments the dual mesh `dualmesh`, the function symbol, and optionally the hodge operator.

### Composing Operators

Decapodes uses the Discrete Exterior Calculus to discretize our differential operators. The DEC is an elegant way of building up more complex differential operators from simpler ones. To demonstrate, we will define the Δ operator by building it up with matrix multiplication of simpler operators. Since most operators in the DEC are matrices, most simulations consist mainly of matrix-vector multiplications, and are thus very fast.

If this code seems too low level, do not worry. Decapodes defines and caches for you many differential operators behind the scenes, so you do not have to worry about defining your own.

```
lap_mat = dec_hodge_star(1,dualmesh) * dec_differential(0,dualmesh) * dec_inv_hodge_star(0,dualmesh) * dec_dual_derivative(0,dualmesh)

function generate(dualmesh, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :Δ => x -> begin
      lap_mat * x
    end
  end
  return (args...) -> op(args...)
end
```
