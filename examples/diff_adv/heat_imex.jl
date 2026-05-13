using CombinatorialSpaces
using ComponentArrays
using Decapodes
using DiagrammaticEquations
using GeometryBasics: Point3
using OrdinaryDiffEq

Point3D = Point3{Float64}

# Heat equation with linear decay: ∂ₜ(C) == κ*Δ(C) + S*C, S < 0.
#   Implicit (stiff) part:     ∂ₜ(C) == κ*Δ(C)
#   Explicit (non-stiff) part: ∂ₜ(C) == S*C
heat_implicit = @decapode begin
  C::Form0
  κ::Constant
  ∂ₜ(C) == κ * Δ(C)
end

heat_explicit = @decapode begin
  C::Form0
  S::Constant
  ∂ₜ(C) == S * C
end

grid_size = 10
s  = triangulated_grid(grid_size, grid_size, 1, 1, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

split_sim = eval_split(heat_implicit, heat_explicit, preallocate=false)
f_implicit, f_explicit = split_sim(sd, nothing, DiagonalHodge())

C_vals = map(sd[:point]) do (x, y)
  sin(pi * x / grid_size) * sin(pi * y / grid_size)
end

u₀ = ComponentArray(C = C_vals)
p  = (κ = 0.1, S = -1.0)

prob = SplitODEProblem(f_implicit, f_explicit, u₀, (0.0, 0.5), p)
soln = solve(prob, KenCarp4())
soln.retcode != :Unstable || error("Solver was not stable")
