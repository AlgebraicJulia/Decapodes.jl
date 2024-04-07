using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using Catlab
using CombinatorialSpaces
using Distributions
using CairoMakie
using JLD2
using LinearAlgebra
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using GeometryBasics: Point2
Point2D = Point2{Float64}

Klausmeier = @decapode begin
  (n,w)::DualForm0
  dX::Form1
  (a,ν,m)::Constant

  ∂ₜ(w) == a - w - w * n^2 + ν * L(dX, w)
  ∂ₜ(n) == w * n^2 - m*n + Δ(n)
end

Klausmeier[9, :type] = :DualForm0
Klausmeier[10, :type] = :DualForm0
Klausmeier[15, :type] = :DualForm0

sim = eval(gensim(Klausmeier, dimension=1))

function circle(n, c)
  s = EmbeddedDeltaSet1D{Bool, Point2D}()
  map(range(0, 2pi - (pi/(2^(n-1))); step=pi/(2^(n-1)))) do t
    add_vertex!(s, point=Point2D(cos(t),sin(t))*(c/2pi))
  end
  add_edges!(s, 1:(nv(s)-1), 2:nv(s))
  add_edge!(s, nv(s), 1)
  sd = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s)
  subdivide_duals!(sd, Circumcenter())
  s,sd
end
s,sd = circle(7, 500)

lap_mat = hodge_star(1,sd) * d(0,sd) * inv_hodge_star(0,sd) * dual_derivative(0,sd)
function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :Δ => x -> begin
      lap_mat * x
    end
    :^ => (x, y) -> x .^ y
  end
  return (args...) -> op(args...)
end

fₘ = sim(sd, generate, DiagonalHodge())

n_dist = Normal(pi)
n = [pdf(n_dist, t)*(√(2pi))*7.2 + 0.08 - 5e-2 for t in range(0,2pi; length=ne(sd))]

w_dist = Normal(pi, 20)
w = [pdf(w_dist, t) for t in range(0,2pi; length=ne(sd))]

dX = sd[:length]

u₀ = ComponentArray(n = n, w = w, dX = dX)

cs_ps = (m = 0.45,
         a = 0.94,
         ν = 182.5)

tₑ = 600.0
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), cs_ps)
@info "Solving"
soln = solve(prob, Tsit5())
                 