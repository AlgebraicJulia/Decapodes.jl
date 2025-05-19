# Parameter Calibration for Halfar's Glacial Flow Model

```@setup INFO
include(joinpath(Base.@__DIR__, "..", "..", "docinfo.jl"))
info = DocInfo.Info()
```

Let's see how to calibrate a glacial flow Decapode model's parameters to fit some data. We want to solve the inverse problem, i.e. given a model and some data, find a set of parameters that causes the output of the model to match the given data as closely as possible. 

We'll be using a 2D version of the Halfar glacial flow model, for more explanation see the [glacial flow docs page](../ice_dynamics/ice_dynamics.md).

In order to set up the inverse problem, we first need our model and some reference data. So, we'll set up the 2D glacial flow model and get some data from it. In this case we'll be fitting our model parameters to data from the model itself. In general the data we want to fit to will not be from the model, but for demonstration purposes this works well. 

``` @example Calibration
using Catlab
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes
using CairoMakie
using ComponentArrays
using GeometryBasics: Point2, Point3
using JLD2
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using Statistics
using BenchmarkTools
using SparseConnectivityTracer
using SparseMatrixColorings
using ADTypes
using NaNMath
using Optimization
using SciMLSensitivity
using Optimization
using OptimizationOptimJL
using Optim
Point2D = Point2{Float64}
Point3D = Point3{Float64}

halfar_eq2 = @decapode begin
    h::Form0
    Γ::Form1
    n::Constant

    ḣ == ∂ₜ(h)
    ḣ == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯(d(h)))^(n - 1)) * avg₀₁(h^(n + 2)))
end

glens_law = @decapode begin
    #Γ::Form0
    Γ::Form1
    (A, ρ, g, n)::Constant

    Γ == (2 / (n + 2)) * A * (ρ * g)^n
end

ice_dynamics_composition_diagram = @relation () begin
    dynamics(Γ, n)
    stress(Γ, n)
end

ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
    [Open(halfar_eq2, [:Γ, :n]),
        Open(glens_law, [:Γ, :n])])

ice_dynamics = apex(ice_dynamics_cospan)

ice_dynamics2D = expand_operators(ice_dynamics)
infer_types!(ice_dynamics2D)
resolve_overloads!(ice_dynamics2D)

function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
        :♯ => begin
            sharp_mat = ♯_mat(sd, AltPPSharp())
            x -> sharp_mat * x
        end
        :mag => x -> norm.(x)
        x => error("Unmatched operator $my_symbol")
    end
    return (args...) -> op(args...)
end
```
Now we define the mesh and the dual mesh that the glacial flow equations will be solved on, along with initial conditions on the mesh and parameters for the model.

```@example Calibration
s2D = triangulated_grid(10_000, 10_000, 500, 500, Point3D)
sd2D = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s2D)
subdivide_duals!(sd2D, Barycenter())

h₀2D = map(point(s2D)) do (x, y)
    (7072 - ((x - 5000)^2 + (y - 5000)^2)^(1 / 2)) / 9e3 + 10
end

u02D = ComponentArray(dynamics_h=h₀2D)

n = 3
ρ = 910
g = 9.8
A = 1e-3

constants_and_parameters = ComponentArray(
    n=n,
    stress_ρ=ρ,
    stress_g=g,
    stress_A=A)
```

Now we have everything we need to generate the function that will be used in the ODE solver. Thinking ahead however, the optimization routines we want to use might try non-physical parameters, e.g. maybe a set of parameters the optimizer tries will cause the ice height solution to become complex. In this case, certain parameter sets might cause an error to occur, e.g. an out of domain error could be thrown if the square root function is called on a negative number. If an error occurs, the optimization process will stop. But we would like it if when a non-physical set of parameters is used, the optimization routine rejects it and keeps going. For that we can use the [`NaNMath` library](https://github.com/JuliaMath/NaNMath.jl). When functions in the `NaNMath` library are called with values that would usually produce an out of domain error, they return a `NaN` value instead. This allows for the optimization routine to continue, while disregarding that set of parameters. 

```@example Calibration

decapode_code =
    quote
        let
            ^(x,y) = NaNMath.pow(x,y)
            sqrt(x) = NaNMath.sqrt(x)
            log(x) = NaNMath.log(x)
            
            $(gensim(ice_dynamics2D, dimension = 2, preallocate = false))
        end
    end

f_eval_2D = eval(decapode_code)

f_2D = f_eval_2D(sd2D,generate)
```

Now we can solve the problem and generate some reference data.

```@example Calibration
tₑ = 8e3

data_prob = ODEProblem{true,SciMLBase.FullSpecialize}(f_2D, u02D, (0, tₑ), constants_and_parameters)
decapode_sol = solve(data_prob, Rodas5P())

reference_dat = last(decapode_sol).dynamics_h
```

In order to optimize the parameters of the model, we compare the output of the model with certain parameters against the reference data. When the output of the model is as close as possible to the reference data, we'll have a set of parameters that when used in the model will reproduce the reference data. The difference between the model output and the reference data can be calculated in many ways, in this case we'll use the sum of squares. 

```@example Calibration
function loss(u) #only compares last time step
    newp = ComponentArray(n=n, stress_ρ=u[1], stress_g=g, stress_A=A)
    prob = remake(data_prob, p=newp)
    sol = solve(prob, Rodas5P(), sensealg=GaussAdjoint())
    current_dat = last(sol).dynamics_h
    sum(abs2, reference_dat .- current_dat)
end
```

## Sparsity Detection and Jacobian Coloring
Before we get to setting up and solving the optimization problem, let's see if there are any easy things we can do to make sure we have good performance. The optimization routine will run the `loss` function with many different sets of parameters, each time checking to see how close the output of the model with that set of parameters is. If we want to speed up the optimization process, one way to do that is to speed up the solution of the ordinary differential equation involved. 

Decapodes used a semi-discretization method to discretize the spatial variables of the model. This means that the system of ordinary differential equations that was generated has a state associated with every vertex of the mesh. We chose our mesh to be 21 by 21, so there are 441 vertices in our mesh, which means our ODE has 441 state variables. Since the system of equations might be [stiff](https://en.wikipedia.org/wiki/Stiff_equation), we want to use an implicit method, which means that for every time step the ODE solver takes, it needs to calculate the Jacobian of `f_2D` with respect to the state variables. Because we have have 441 inputs, and 441 outputs, the Jacobian matrix has a total of 194_481 entries. This seems like a lot, but actually, since not every output will depend on every input, the Jacobian will have many entries that will always be zero. In other words the Jacobian will be sparse. This means that only a few entries actually have to be calculated. We can use the package `SparseConnectivityTracer.jl` to find the sparsity pattern of the Jacobian.  

```@example Calibration
jac_sparsity_2D = jacobian_sparsity((du,u) -> f_2D(du,u,constants_and_parameters,0.0),u02D, u02D, TracerLocalSparsityDetector())
```

We can see that only 7597 entries are actually ever non-zero. If we use a sparsity aware automatic differentiation system to calculate the Jacobian while solving the ODE, we can cut down the number of calculations needed by two orders of magnitude.  

There's another way we can speed up the calculation of the Jacobian. Forward mode automatic differentiation systems use the Jacobian vector product (JVP) as their base operation. Meaning that the Jacobian is built column by column using many JVP operations. One way to cut down on the number of JVP operations needed to get the full Jacobian is to use Jacobian coloring. In many cases, the sparse Jacobian will have many columns that do not share any non-zero entries. In that case the columns can be "compressed" in to one column. Now all of the original columns of the Jacobian that were compressed in to one column can be calculated using one JVP operation. This can dramatically cut down the number of operations needed to evaluate a Jacobian. It's called a Jacobian coloring because seperating a matrix in to columns that do not share a non-zero entry is equivalent to a graph coloring problem. Here we can see that using a greedy coloring algorithm from `SparseMatrixColorings.jl`, the Jacobian for our function only has 26 colors.  

```@example Calibration
col_prob = ColoringProblem()
algo = GreedyColoringAlgorithm()

jac_colors_2D = coloring(jac_sparsity_2D,col_prob,algo)

ncolors(jac_colors_2D)
```

Now we can take the coloring pattern we found and see what it looks like when it's applied to the sparsity pattern. 

```@example Calibration
compress(jac_sparsity_2D,jac_colors_2D)

sparse(compress(jac_sparsity_2D,jac_colors_2D)) # hide
```
We can see that this matrix has only 26 columns. This means that when the ODE solver calculates the derivative of the function, it can do it in 26 JVP operations instead of 441, a huge difference! Once the Jacobian is calculated in this compressed form, it can then be "decompressed" using the original coloring pattern. 

Now we can look at how much exactly this has sped up the solving process.

First, let's solve the ODE with sparsity not used.

```@example Calibration
no_sparse_prob_2D = ODEProblem(f_2D, u02D, (0, tₑ), constants_and_parameters)
no_sparse_soln_2D, exec_time_seconds, _, _, _ = @btimed solve(no_sparse_prob_2D, Rodas5P(autodiff = AutoForwardDiff()))
no_sparse_soln_2D.retcode, exec_time_seconds
```

Now the same problem but with the sparsity pattern and Jacobian coloring taken in to account. 
```@example Calibration
sparse_f_2D = ODEFunction(f_2D, sparsity = jac_sparsity_2D, colorvec = column_colors(jac_colors_2D))
sparse_prob_2D = ODEProblem(sparse_f_2D,u02D,(0,tₑ), constants_and_parameters)
sparse_soln_2D, exec_time_seconds, _, _, _ = @btimed solve(sparse_prob_2D, Rodas5P(autodiff = AutoForwardDiff()))
sparse_soln_2D.retcode, exec_time_seconds
```

We can see quite a performance increase!

For more details on sparsity detection and Jacobian coloring, see the preprint [Sparser, Better, Faster, Stronger: Efficient Automatic Differentiation for Sparse Jacobians and Hessians](https://arxiv.org/abs/2501.17737) from A. Hill and G. Dalle.

## Calibration

Now that we know we can take advantage of sparsity detection and Jacobian coloring, we actually have to tweak our loss function to use them.

```@example Calibration
function loss(u) #only compares last time step
    newp = ComponentArray(n=n, stress_ρ=u[1], stress_g=g, stress_A=A)
    prob = remake(sparse_prob_2D, p=newp)
    sol = solve(prob, Rodas5P())
    current_dat = last(sol).dynamics_h
    sum(abs2, reference_dat .- current_dat)
end
```

And now to calibrate the $\rho$ parameter, all we can use `Optimization.jl` to set up an `OptimizationFunction` with our loss function, and an `OptimizationProblem` where the initial guess is `100.0`. The specific optimization routine used is the [`Optim.jl` LFBGS](https://julianlsolvers.github.io/Optim.jl/stable/algo/lbfgs/) routine with default settings. This is a gradient based quasi-Newton method that uses automatic differentiation to find the gradient, and uses it to estimate the Hessian. Forward mode automatic differentiation from [`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/stable/) is used to differentiate the loss function to provide the gradients used by LFBGS. Keep in mind that the value of $\rho$ used to generate the reference data was 910. If the optimization works we would expect to get a value very close to that.

```@example Calibration
optf = OptimizationFunction((x, p) -> loss(x), AutoForwardDiff())
optprob = Optimization.OptimizationProblem(optf, [100.0])
optsol = Optimization.solve(optprob, LBFGS())
```

As we can see the optimization routine was able to find the parameter that produced the reference data quite closely. 

