using Decapodes
using Catlab
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
Point2D = Point2{Float64}
Point3D = Point3{Float64}

using CUDA
using CUDA.CUSPARSE

rect = triangulated_grid(100, 100, 0.2, 0.2, Point3D);
d_rect = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(rect);
subdivide_duals!(d_rect, Circumcenter());

CahnHillard = @decapode begin
    C::Form0
    ∂ₜ(C) == 0.5 * Δ(C.^3 - C - 0.5 * Δ(C))
end

sim = eval(gensim(CahnHillard, code_target=CUDATarget()))
  
fₘ = sim(d_rect, nothing)

C = CUDA.rand(Float64, nv(d_rect)) * 2 .- 1

u₀ = ComponentArray(C=C)

tₑ = 60

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4))
soln = solve(prob, Tsit5());
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ))
soln = solve(prob, Tsit5(), saveat=0.1, progress=true);
@info("Done")

begin
    frames = 200
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    msh = CairoMakie.mesh!(ax, rect, color=Array(soln(0).C), colormap=:jet, colorrange=extrema(Array(soln(0).C)))
    Colorbar(fig[1,2], msh)
    CairoMakie.record(fig, "CahnHillard.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
      msh.color = Array(soln(t).C)
    end
  end