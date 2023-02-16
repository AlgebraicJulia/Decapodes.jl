using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using GLMakie
using Logging
using JLD2

using GeometryBasics: Point2
Point2D = Point2{Float64}

Brusselator = SummationDecapode(parse_decapode(
quote
  # Values living on vertices.
  (U, V, F)::Form0{X} # State variables.
  (U2V, One)::Form0{X} # Named intermediate variables.
  (U̇, V̇)::Form0{X} # Tangent variables.
  # Scalars.
  (fourfour, threefour, α)::Constant{X}
  # A named intermediate variable.
  U2V == (U .* U) .* V
  # Specify how to compute the tangent variables.
  U̇ == One + U2V - (fourfour * U) + (α * Δ(U)) + F
  V̇ == (threefour * U) - U2V + (α * Δ(U))
  # Associate tangent variables with a state variable.
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end))
# Visualize. You must have graphviz installed.
to_graphviz(Brusselator)

# We resolve types of intermediate variables using sets of rules.
bespoke_op1_inf_rules = [
  (src_type = :Form0, tgt_type = :infer, replacement_type = :Form0, op = :Δ)]

bespoke_op2_inf_rules = [
  (proj1_type = :Form0, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :.*),
  (proj1_type = :Form0, proj2_type = :Parameter, res_type = :infer, replacement_type = :Form0, op = :*),
  (proj1_type = :Parameter, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :*)]

infer_types!(Brusselator,
    vcat(bespoke_op1_inf_rules, op1_inf_rules_2D),
    vcat(bespoke_op2_inf_rules, op2_inf_rules_2D))
# Visualize. Note that variables now all have types.
to_graphviz(Brusselator)

# Resolve overloads. i.e. ~dispatch
resolve_overloads!(Brusselator)
# Visualize. Note that functions are renamed.
to_graphviz(Brusselator)

# TODO: Create square domain of approximately 32x32 vertices.
s = loadmesh(Rectangle_30x10())
scaling_mat = Diagonal([1/maximum(x->x[1], s[:point]),
                        1/maximum(x->x[2], s[:point]),
                        1.0])
s[:point] = map(x -> scaling_mat*x, s[:point])
orient!(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())

# Define how operations map to Julia functions.
hodge = GeometricHodge()
Δ₀ = δ(1, sd, hodge=hodge) * d(0, sd)
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    # The Laplacian operator on 0-Forms is the codifferential of
    # the exterior derivative. i.e. dδ
    :Δ₀ => x -> Δ₀ * x
    :.* => (x,y) -> x .* y
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

# Create initial data.
@assert all(map(sd[:point]) do (x,y)
  0.0 ≤ x ≤ 1.0 && 0.0 ≤ y ≤ 1.0
end)

U = map(sd[:point]) do (_,y)
  22 * (y *(1-y))^(3/2)
end

V = map(sd[:point]) do (x,_)
  27 * (x *(1-x))^(3/2)
end

F₁ = map(sd[:point]) do (x,y)
 (x-0.3)^2 + (y-0.6)^2 ≤ (0.1)^2 ? 5.0 : 0.0
end

F₂ = zeros(nv(sd))

One = ones(nv(sd))

constants = (
  fourfour = 4.4,
  threefour = 3.4,
  α = 0.001)

# Generate the simulation.
gensim(expand_operators(Brusselator))
sim = eval(gensim(expand_operators(Brusselator)))
fₘ = sim(sd, generate)

# Create problem and run sim for t ∈ [0,1.1).
# Map symbols to data.
u₀ = construct(PhysicsState, [VectorForm(U), VectorForm(V), VectorForm(F₁), VectorForm(One)], Float64[], [:U, :V, :F, :One])

# Visualize the initial conditions.
# If GLMakie throws errors, then update your graphics drivers,
# or use an alternative Makie backend like CairoMakie.
fig_ic = GLMakie.Figure()
p1 = GLMakie.mesh(fig_ic[1,2], s, color=findnode(u₀, :U), colormap=:jet)
p2 = GLMakie.mesh(fig_ic[1,3], s, color=findnode(u₀, :V), colormap=:jet)
p3 = GLMakie.mesh(fig_ic[1,4], s, color=findnode(u₀, :F), colormap=:jet)

tₘ = 1.1

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₘ), constants)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator_middle.jld2" soln

GLMakie.mesh(s, color=findnode(soln(tₘ), :U), colormap=:plasma)

begin # BEGIN Gif creation
times = range(0.0, tₘ, length=75)
colors_U = [findnode(soln(t), :U) for t in times]
colors_V = [findnode(soln(t), :V) for t in times]
# Initial frame
fig = GLMakie.Figure(resolution = (1200, 800))
p1 = GLMakie.mesh(fig[1,2], s, color=colors_U[1], colormap=:jet, colorrange=extrema(colors_U[1]))
p2 = GLMakie.mesh(fig[1,4], s, color=colors_V[1], colormap=:jet, colorrange=extrema(colors_V[1]))
ax1 = Axis(fig[1,2], width = 400, height = 400)
ax2 = Axis(fig[1,4], width = 400, height = 400)
hidedecorations!(ax1)
hidedecorations!(ax2)
hidespines!(ax1)
hidespines!(ax2)
Colorbar(fig[1,1])
Colorbar(fig[1,5])
Label(fig[1,2,Top()], "U")
Label(fig[1,4,Top()], "V")
lab1 = Label(fig[1,3], "")

# Animation
using Printf
record(fig, "brusselator_middle.gif", range(0.0, tₘ; length=75); framerate = 30) do t
    p1.plot.color = findnode(soln(t), :U)
    p2.plot.color = findnode(soln(t), :V)
    lab1.text = @sprintf("%.2f",t)
end

end # END Gif creation

# Create problem and run sim for t ∈ [1.1,11.5].
# Map symbols to data.
u₁ = construct(PhysicsState,
  [findnode(soln(tₘ), :U),
   findnode(soln(tₘ), :V),
   VectorForm(F₂),
   VectorForm(One)], Float64[], [:U, :V, :F, :One])

# Visualize the initial conditions.
# If GLMakie throws errors, then update your graphics drivers,
# or use an alternative Makie backend like CairoMakie.
fig_ic = GLMakie.Figure()
p1 = GLMakie.mesh(fig_ic[1,2], s, color=findnode(u₁, :U), colormap=:jet)
p2 = GLMakie.mesh(fig_ic[1,3], s, color=findnode(u₁, :V), colormap=:jet)
p3 = GLMakie.mesh(fig_ic[1,4], s, color=findnode(u₁, :F), colormap=:jet)

tₑ = 11.5

@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₁, (0, 1e-4), constants)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ, u₁, (tₘ, tₑ), constants)
soln = solve(prob, Tsit5())
@info("Done")

@save "brusselator_end.jld2" soln

GLMakie.mesh(s, color=findnode(soln(tₘ), :U), colormap=:plasma)

begin # BEGIN Gif creation
num_frames = Int(floor(75 * 10.4 / 1.1))
times = range(tₘ, tₑ, length=num_frames)
colors_U = [findnode(soln(t), :U) for t in times]
colors_V = [findnode(soln(t), :V) for t in times]
# Initial frame
fig = GLMakie.Figure(resolution = (1200, 800))
p1 = GLMakie.mesh(fig[1,2], s, color=colors_U[1], colormap=:jet, colorrange=extrema(colors_U[1]))
p2 = GLMakie.mesh(fig[1,4], s, color=colors_V[1], colormap=:jet, colorrange=extrema(colors_V[1]))
ax1 = Axis(fig[1,2], width = 400, height = 400)
ax2 = Axis(fig[1,4], width = 400, height = 400)
hidedecorations!(ax1)
hidedecorations!(ax2)
hidespines!(ax1)
hidespines!(ax2)
Colorbar(fig[1,1])
Colorbar(fig[1,5])
Label(fig[1,2,Top()], "U")
Label(fig[1,4,Top()], "V")
lab1 = Label(fig[1,3], "")

# Animation
using Printf
record(fig, "brusselator_end.gif", range(tₘ, tₑ; length=num_frames); framerate = 30) do t
    p1.plot.color = findnode(soln(t), :U)
    p2.plot.color = findnode(soln(t), :V)
    lab1.text = @sprintf("%.2f",t)
end

end # END Gif creation
