# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using DiagrammaticEquations.Deca

# External Dependencies
using ComponentArrays
using MLStyle
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using SparseArrays
using Statistics
using CairoMakie
using CUDA
using CUDA.CUSPARSE
using BenchmarkTools
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ  * d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end

# Equation 1 from Glen, J. W. THE FLOW LAW OF ICE: A discussion of the
# assumptions made in glacier theory, their experimental foundations and
# consequences. (1958)
glens_law = @decapode begin
  Γ::Form1
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end

#####################
# Compose the model #
#####################

ice_dynamics_composition_diagram = @relation () begin
  dynamics(Γ,n)
  stress(Γ,n)
end

# Plug in our Decapodes to the composition pattern.
ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:Γ,:n]),
   Open(glens_law, [:Γ,:n])])

ice_dynamics2D = apex(ice_dynamics_cospan)

# Interpret this multiphysics diagram in the 2D exterior calculus.

s′ = triangulated_grid(60_000,100_000,2_000,2_000,Point3D)
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())
x̄ = mean(p -> p[1], point(s))
ȳ = mean(p -> p[2], point(s))

# These are the initial conditions to the Halfar Dome test case that the
# Community Ice Sheet Model uses.
R₀ = 60_000 * sqrt(0.125)
H = 2_000 * sqrt(0.125)

n = 3
g = 9.8101
ρ = 910
alpha = 1/9
beta = 1/18
flwa = 1e-16
A = fill(1e-16, ne(s))

Gamma = 2.0/(n+2) * flwa * (ρ * g)^n
t0 = (beta/Gamma) * (7.0/4.0)^3 * (R₀^4 / H^7)

# This is the analytic solution for comparison.
# It is ported over from the CISM code for comparison's sake,
# and we will use it to set initial conditions.
function height_at_p(x,y,t)
  tr = (t + t0) / t0
  r = sqrt((x - x̄)^2 + (y - ȳ)^2)
  r = r/R₀
  inside = max(0.0, 1.0 - (r / tr^beta)^((n+1.0) / n))
  H * inside^(n / (2*n + 1)) / tr^alpha
end

# Set the initial conditions for ice sheet height:
# Ice height is a primal 0-form. i.e. valued at vertices.
h₀ = map(x -> height_at_p(x[1], x[2], 0), point(s′))
fig = mesh(s′, color=h₀, colormap=:jet)

# Store these values to be passed to the solver.
u₀ = ComponentArray(dynamics_h = CuVector{Float64}(h₀))
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = CuVector{Float64}(A))

#############################################
# Define how symbols map to Julia functions #
#############################################

# This sharp operator, ♯, is scheduled to be upstreamed.
function generate(sd, my_symbol; hodge=GeometricHodge())
  # We pre-allocate matrices that encode differential operators.
  op = @match my_symbol begin
    :♯ => begin 
      # TODO: For some reason this works as a dense CuArray but not when sparse
      ♯_m = CuArray(♯_mat(sd, LLSDDSharp()))
      x -> ♯_m * x
    end
    :mag => x -> begin
      CUDA.norm.(x)
    end
    :^ => (x,y) -> begin
      x .^ y
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

#######################
# Generate simulation #
#######################

sim = eval(gensim(ice_dynamics2D, code_target=gen_CUDA()))
fₘ = sim(s, generate)

# Pre-compile simulation

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
# We run for a short timespan to pre-compile.
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

# Run simulation
tₑ = 200

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")


# Plot the final conditions
function plot_final_conditions()
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1],
    title="Modeled thickness (m) at time 200.0",
    aspect=0.6)
  msh = mesh!(ax, s′, color=Array(soln(200.0).dynamics_h), colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
fig = plot_final_conditions()
save("ice_numeric_solution.png", fig)

# Plot the final conditions according to the analytic solution.
function plot_analytic()
  hₐ = map(x -> height_at_p(x[1], x[2], 200.0), point(s′))
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1],
    title="Analytic thickness (m) at time 200.0",
    aspect=0.6)
  msh = mesh!(ax, s′, color=hₐ, colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
fig = plot_analytic()
save("ice_analytic_solution.png", fig)

# Plot the error.
function plot_error()
  hₐ = map(x -> height_at_p(x[1], x[2], 200.0), point(s′))
  h_diff = Array(soln(tₑ).dynamics_h) - hₐ
  extrema(h_diff)
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1],
    title="Modeled thickness - Analytic thickness at time 200.0",
    aspect=0.6)
  msh = mesh!(ax, s′, color=h_diff, colormap=:jet)
  Colorbar(fig[1,2], msh)
  fig
end
fig = plot_error()
save("ice_error.png", fig)

# Compute max absolute error:
hₐ = map(x -> height_at_p(x[1], x[2], 200.0), point(s′))
h_diff = Array(soln(tₑ).dynamics_h) - hₐ
maximum(abs.(h_diff))

# Compute RMSE not considering the "outside".
hₐ = map(x -> height_at_p(x[1], x[2], 200.0), point(s′))
nonzeros = findall(!=(0), hₐ)
h_diff = Array(soln(tₑ).dynamics_h) - hₐ
rmse = sqrt(sum(map(x -> x*x, h_diff[nonzeros])) / length(nonzeros))

# Compute RMSE of the entire domain.
hₐ = map(x -> height_at_p(x[1], x[2], 200.0), point(s′))
h_diff = Array(soln(tₑ).dynamics_h) - hₐ
rmse = sqrt(sum(map(x -> x*x, h_diff)) / length(h_diff))

# Create a gif
begin
  frames = 300
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s′, color=Array(soln(0).dynamics_h), colormap=:jet, colorrange=extrema(Array(soln(0).dynamics_h)))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "Shallow_Ice_GPU.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
      msh.color = Array(soln(t).dynamics_h)
  end
end