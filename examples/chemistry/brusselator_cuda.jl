using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using CairoMakie
using GeometryBasics: Point2, Point3
using CUDA
using CUDA.CUSPARSE
Point2D = Point2{Float64}
Point3D = Point3{Float64}

Brusselator = @decapode begin
  # Values living on vertices.
  (U, V)::Form0{X} # State variables.
  (U2V)::Form0{X} # Named intermediate variables.
  (U̇, V̇)::Form0{X} # Tangent variables.
  # Scalars.
  (α)::Constant{X}
  F::Parameter{X}
  # A named intermediate variable.
  U2V == (U .* U) .* V
  # Specify how to compute the tangent variables.
  U̇ == 1 + U2V - (4.4 * U) + (α * Δ(U)) + F
  V̇ == (3.4 * U) - U2V + (α * Δ(U))
  # Associate tangent variables with a state variable.
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

gensim(Brusselator, code_target = gen_CUDA())

# TODO: Create square domain of approximately 32x32 vertices.
# s = loadmesh(Rectangle_30x10())
# Visualize the mesh.
s = triangulated_grid(200, 200, 5, 5, Point3D)
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0])
s[:point] = map(x -> scaling_mat*x, s[:point])
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

U = CuVector{Float64}(map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end)

V = CuVector{Float64}(map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end)

F₁ = CuVector{Float64}(map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end)

F₂ = CUDA.zeros(Float64, nv(sd))

constants_and_parameters = (
  α = 0.001,
  F = t -> t ≥ 1.1 ? F₂ : F₁)

sim = eval(gensim(Brusselator, code_target = gen_CUDA()))
fₘ = sim(sd, generate, DiagonalHodge())

# Create problem and run sim for t ∈ [0,tₑ).
# Map symbols to data.
u₀ = ComponentArray(U=U, V=V)

tₑ = 25

begin
  @info("Precompiling Solver")
  prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
  soln = solve(prob, Tsit5())
  soln.retcode != :Unstable || error("Solver was not stable")
  @info("Solving")
  prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
  soln = solve(prob, Tsit5())
  @info("Done")
end


#=
begin
    frames = 300
    fig = Figure()
    ax = CairoMakie.Axis(fig[1,1])
    msh = CairoMakie.mesh!(ax, s, color=soln(0).U, colormap=:jet, colorrange=extrema(soln(0).U))
    Colorbar(fig[1,2], msh)
    CairoMakie.record(fig, "Brusselator_GPU.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
        msh.color = soln(t).U
    end
end
=#