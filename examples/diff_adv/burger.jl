# AlgebraicJulia Dependencies
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using Catlab
using CombinatorialSpaces

# External Dependencies
using Logging: global_logger
## using TerminalLoggers: TerminalLogger
## global_logger(TerminalLogger())
using GeometryBasics: Point2
Point2D = Point2{Float64}
using Distributions
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using CairoMakie
import CairoMakie: Axis
using ComponentArrays

# Represent component Decapodes.
Diffusion = @decapode begin
  C::Form0
  ϕ::Form1
  ν::Constant

  ## Fick's first law
  ϕ == ν * d(C)
end

Advection = @decapode begin
  C::Form0
  (V, ϕ)::Form1

  ϕ == ∧₀₁(C,V)
end

Lie = @decapode begin
  C::Form0
  V::Form1
  dX::Form1

  V == ∘(⋆,⋆)(C ∧ dX)
end

Superposition = @decapode begin
  (C, Ċ)::Form0
  (ϕ, ϕ₁, ϕ₂)::Form1

  ϕ == ϕ₁ + ϕ₂
  Ċ == ∘(⋆,d,⋆)(ϕ)
  ∂ₜ(C) == Ċ
end

# Compose physics.
compose_burger = @relation () begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  lie(C, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

to_graphviz(compose_burger, box_labels=:name, junction_labels=:variable, prog="circo")

Burger_cospan = oapply(compose_burger,
               [Open(Diffusion,     [:C, :ϕ]),
                Open(Advection,     [:C, :ϕ, :V]),
                Open(Lie,           [:C, :V]),
                Open(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])])
Burger = apex(Burger_cospan)

# Specify semantics of the 1D DEC.
# i.e. Declare these dynamics are happening on a line.
Burger = expand_operators(Burger)
infer_types!(Burger, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(Burger, op1_res_rules_1D, op2_res_rules_1D)

to_graphviz(Burger)

# Create mesh.
# This is a line. This could be a helper function.
s = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s, 1000, point=Point2D.(1:1000,0))
add_edges!(s, 1:(nv(s)-1), 2:nv(s))
sd = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s)
subdivide_duals!(sd, Circumcenter())

# Set initial conditions and constants.
c_dist = MvNormal([500, 5], [10.5, 10.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in point(sd)]
dX = ones(ne(sd))

u₀ = ComponentArray(C=c, lie_dX=dX)

cs_ps = (diffusion_ν = 0.0005,)

# Describe mappings from symbols to discrete differential operators.
function generate(sd, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    ## Specify which wedge product to use.
    ## This should probably be the default.
    :∧₀₁ => (x,y) -> begin
      ∧(Tuple{0,1},sd,x,y)
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

# Generate simulation.
sim = eval(gensim(Burger, dimension=1))
fₘ = sim(sd, generate, DiagonalHodge())

# Run simulation.
tₑ = 1e5
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), cs_ps)
sol = solve(prob, Tsit5(), progress=true, progress_steps=1)

# Visualize initial and final conditions.
lines(map(x -> x[1], point(sd)), sol(0.0).C)
lines!(map(x -> x[1], point(sd)), sol(tₑ).C)

# Animate the dynamics.
times = range(0.0, tₑ, length=150)
colors = [sol(t).C for t in times]

frames = 100
fig = Figure(resolution = (800, 800))
ax1 = Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(sd))))
ylims!(ax1, extrema(sol(0.0).C))
Label(fig[1,1,Top()], "Speed C")
Label(fig[2,1,Top()], "Line plot of speed of fluid along the linear domain, every $(tₑ/frames) time units")

record(fig, "burger_low_diff.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(sd)), sol(t).C)
end
