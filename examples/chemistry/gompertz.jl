using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using Distributions
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using CairoMakie
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

function show_heatmap(Cdata)
  heatmap(reshape(Cdata, (floor(Int64, sqrt(length(Cdata))), floor(Int64, sqrt(length(Cdata))))))
end

s = triangulated_grid(50,50,0.2,0.2,Point2D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point2D}(s)
subdivide_duals!(sd, Circumcenter())

constants_and_parameters = (
  Dif = 0.005,
  Kd = 0.5,
  Cmax = 10
)

# Eq. 34 from Yi et al.
# A Review of Mathematical Models for Tumor Dynamics and Treatment Resistance
# Evolution of Solid Tumors,
# with f given as logistic growth. (Eq. 5)
Logistic = @decapode begin
  C::Form0
  (Dif, Kd, Cmax)::Constant

  fC == C * (1 - C / Cmax)
  ∂ₜ(C) == Dif * Δ(C) + fC - Kd * C 
end

sim = eval(gensim(Logistic))

fₘ = sim(sd, nothing, DiagonalHodge())

# "The model ... considers an equivalent radially symmetric tumour"
# - Murray J.D., Glioblastoma brain tumours
c_dist  = MvNormal([25, 25], 2)
C = 100 * [pdf(c_dist, [p[1], p[2]]) for p in sd[:point]]

u₀ = ComponentArray(C=C)
tₑ = 15.0

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

show_heatmap(soln(0).C)
show_heatmap(soln(tₑ).C)

Gompertz = @decapode begin
  C::Form0
  (Dif, Kd, Cmax)::Constant

  fC == C * ln(Cmax / C)
  ∂ₜ(C) == Dif * Δ(C) + fC - Kd * C 
end

sim = eval(gensim(Gompertz))

function generate(sd, symbol, hodge=DiagonalHodge())
  op = @match symbol begin
    :ln => (x -> log.(x))
    _ => error("Unmatched operator $my_symbol")
  end
  return op
end

fₘ = sim(sd, generate, DiagonalHodge())

constants_and_parameters = (
  Dif = 0.005,
  Kd = 0.5,
  Cmax = 10
)

# "The model ... considers an equivalent radially symmetric tumour"
# - Murray J.D., Glioblastoma brain tumours
c_dist  = MvNormal([25, 25], 2)
C = 100 * [pdf(c_dist, [p[1], p[2]]) for p in sd[:point]]

u₀ = ComponentArray(C=C)
tₑ = 15.0

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

show_heatmap(soln(0).C)
show_heatmap(soln(tₑ).C)