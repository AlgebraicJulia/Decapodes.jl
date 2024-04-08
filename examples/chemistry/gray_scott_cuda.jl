using Catlab
using Catlab.Graphics
using CUDA
using CUDA.CUSPARSE
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

# We use the model equations as stated here and use the initial conditions for
# f, k, rᵤ, rᵥ as listed for experiment 4.
# https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/
GrayScott = @decapode begin
  (U, V)::Form0
  (UV2)::Form0
  (U̇, V̇)::Form0
  (f, k, rᵤ, rᵥ)::Constant
  UV2 == (U .* (V .* V))
  U̇ == rᵤ * Δ(U) - UV2 + f * (1 .- U)
  V̇ == rᵥ * Δ(V) + UV2 - (f + k) .* V
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

s = loadmesh(Rectangle_30x10())
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0])
s[:point] = map(x -> scaling_mat*x, s[:point])
s[:edge_orientation] = false
orient!(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
function generate(sd, my_symbol; hodge=GeometricHodge()) end

U = CuVector{Float64}(map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end)

V = CuVector{Float64}(map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end)

constants_and_parameters = (
  f = 0.024,
  k = 0.055,
  rᵤ = 0.01,
  rᵥ = 0.005)

# Generate the simulation.
sim = eval(gensim(GrayScott, code_target=gen_CUDA()))
fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U,V=V)

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")