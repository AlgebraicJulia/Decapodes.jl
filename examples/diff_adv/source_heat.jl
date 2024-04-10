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

s = triangulated_grid(50,50,0.2,0.2,Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())


function generate(sd, symbol, hodge=DiagonalHodge())
  op = @match symbol begin
    _ => error("Unmatched operator $my_symbol")
  end
  return op
end

Heat = @decapode begin
  C::Form0
  c::Constant
  S::Parameter
  ∂ₜ(C) == 3 * c * Δ(C) + c * S # + L(dX, C)
end

sim = eval(gensim(Heat))

fₘ = sim(sd, generate, DiagonalHodge())

t_dist = MvNormal([25, 25], 1)
S = 50 * [pdf(t_dist, [p[1], p[2]]) for p in sd[:point]]

constants_and_parameters = (
  c = 7,
  S = t -> S * (t % 2 - 0.75)
)

c_dist  = MvNormal([25, 25], 5)
C = 100 * [pdf(c_dist, [p[1], p[2]]) for p in sd[:point]]

u₀ = ComponentArray(C=C)
tₑ = 5

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

show_heatmap(soln(0).C)
show_heatmap(soln(tₑ).C)

begin
  frames = 150
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s, color=soln(0).C, colormap=:jet, colorrange=extrema(soln(0).C))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "Ocillating Heat.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
      msh.color = soln(t).C
  end
end
