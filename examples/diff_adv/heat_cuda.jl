using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using CUDA
# using CairoMakie
using LinearAlgebra
using SparseArrays
Point2D = Point2{Float64}
Point3D = Point3{Float64}
CUDA.allowscalar(false)

rect = loadmesh(Rectangle_30x10())
# rect = triangulated_grid(400, 400, 1, 1, Point3D)
d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(rect)
subdivide_duals!(d_rect, Circumcenter())

Heat = @decapode begin
    U::Form0
    ∂ₜ(U) == 100 * Δ(U)
end

sim = eval(gensim(Heat, code_target=CUDATarget()))

function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      x => error("Unmatched operator $my_symbol")
    end
    return op
  end
  
fₘ = sim(d_rect, generate)

U = CuArray(map(d_rect[:point]) do (x,_)
        return x
    end)

u₀ = ComponentArray(U=U)

constants_and_parameters = ()

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4))
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

CUDA.@time solve(prob, Tsit5());

begin
  frames = 100
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, rect, color=Vector{Float64}(soln(0).U), colormap=:jet, colorrange=extrema(soln(0).U))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "Heat_GPU.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
    msh.color = Vector{Float64}(soln(t).U)
  end
end