println("Loading dependencies")
# AlgebraicJulia Dependencies
using ACSets
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations

# External Dependencies
using CairoMakie
using ComponentArrays
using GeometryBasics
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using Random
using CUDA, CUDA.CUSPARSE
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

Point3D = Point3{Float64};

println("Generating physics")
CahnHilliard = @decapode begin
  C::Form0
  (D, γ)::Constant
  ∂ₜ(C) == D * Δ(C.^3 - C - γ * Δ(C))
end

println("Generating mesh")
s = triangulated_grid(300, 300, 0.2, 0.2, Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool, Float32, Point2{Float32}}(s);
subdivide_duals!(sd, Circumcenter());

println("Mesh elements:") 
println(" Vertices: $(nv(sd))")
println(" Edges: $(ne(sd))")
println(" Triangles: $(ntriangles(sd))")

println(" Dual vertices: $(nparts(sd, :DualV))")
println(" Dual edges: $(nparts(sd, :DualE))")
println(" Dual triangles: $(nparts(sd, :DualTri))")

println("Generating initial conditions")
Random.seed!(0)
C = CUDA.rand(Float32, nv(sd)) * 2 .- 1
u₀ = ComponentArray(C=C)
constants = (D = 0.5f0, γ = 0.5f0);

println("Creating simulation code")
sim = eval(gensim(CahnHilliard, code_target=CUDATarget(), stateeltype=Float32))
fₘ = sim(sd, nothing, DiagonalHodge());

println("Setting up problem")
tₑ = 200
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants)
println("Solving")
soln = solve(prob, Tsit5(), saveat=0.1, progress=true, progress_steps=10_000);
println("Solved")

println("Generating video")
get_data(t) = Array(soln(t).C)

framerate = 30

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=get_data(0), colormap=:jet, colorrange=extrema(get_data(0)))
Colorbar(fig[1,2], msh)

timestamps = range(0, tₑ, length = 5 * framerate)

CairoMakie.record(fig, "CH_CU32.mp4", timestamps) do t
  msh.color = get_data(t)
end
