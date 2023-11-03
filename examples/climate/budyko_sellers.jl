#######################
# Import Dependencies #
#######################

# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes

# External Dependencies
using MLStyle
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
# Uncomment to load GLMakie if your system supports it.
# Otherwise, do using CairoMakie
#using GLMakie
using GeometryBasics: Point2
Point2D = Point2{Float64}

####################
# Define the model #
####################

# ϕ := Latitude
# Tₛ(ϕ,t) := Surface temperature
# Q(ϕ,t) := Insolation
# C(ϕ) := Effective heat capacity
# α(ϕ) := Albedo
# A := Longwave emissions at 0°C
# B := Increase in emissions per degree
# D := Horizontal diffusivity

energy_balance = @decapode begin
  (Tₛ, ASR, OLR, HT)::Form0
  (C)::Constant

  Tₛ̇ == ∂ₜ(Tₛ) 

  Tₛ̇ == (ASR - OLR + HT) ./ C
end
to_graphviz(energy_balance)

absorbed_shortwave_radiation = @decapode begin
  (Q, ASR)::Form0
  α::Constant

  ASR == (1 .- α) .* Q
end
to_graphviz(absorbed_shortwave_radiation)

outgoing_longwave_radiation = @decapode begin
  (Tₛ, OLR)::Form0
  (A,B)::Constant

  OLR == A .+ (B .* Tₛ)
end
to_graphviz(outgoing_longwave_radiation)

heat_transfer = @decapode begin
  (HT, Tₛ)::Form0
  (D,cosϕᵖ,cosϕᵈ)::Constant

  HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
end
to_graphviz(heat_transfer)

insolation = @decapode begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end
to_graphviz(insolation)

to_graphviz(oplus([energy_balance, absorbed_shortwave_radiation, outgoing_longwave_radiation, heat_transfer, insolation]), directed=false)

budyko_sellers_composition_diagram = @relation () begin
  energy(Tₛ, ASR, OLR, HT)
  absorbed_radiation(Q, ASR)
  outgoing_radiation(Tₛ, OLR)
  diffusion(Tₛ, HT, cosϕᵖ)
  insolation(Q, cosϕᵖ)
end
to_graphviz(budyko_sellers_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")

budyko_sellers_cospan = oapply(budyko_sellers_composition_diagram,
  [Open(energy_balance, [:Tₛ, :ASR, :OLR, :HT]),
   Open(absorbed_shortwave_radiation, [:Q, :ASR]),
   Open(outgoing_longwave_radiation, [:Tₛ, :OLR]),
   Open(heat_transfer, [:Tₛ, :HT, :cosϕᵖ]),
   Open(insolation, [:Q, :cosϕᵖ])])

budyko_sellers = apex(budyko_sellers_cospan)
to_graphviz(budyko_sellers)

infer_types!(budyko_sellers, op1_inf_rules_1D, op2_inf_rules_1D)
to_graphviz(budyko_sellers)

resolve_overloads!(budyko_sellers, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(budyko_sellers)

###############################
# Demonstrate storing as JSON #
###############################

write_json_acset(budyko_sellers, "budyko_sellers.json")
# When reading back in, we specify that all attributes are "Symbol"s.
budyko_sellers2 = read_json_acset(SummationDecapode{Symbol,Symbol,Symbol}, "budyko_sellers.json")
# Or, you could choose to interpret the data as "String"s.
budyko_sellers3 = read_json_acset(SummationDecapode{String,String,String}, "budyko_sellers.json")

###################
# Define the mesh #
###################

s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s′, 30, point=Point2D.(range(-π/2 + π/32, π/2 - π/32, length=30), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())

########################################################
# Define constants, parameters, and initial conditions #
########################################################

# This is a primal 0-form, with values at vertices.
cosϕᵖ = map(x -> cos(x[1]), point(s′))
# This is a dual 0-form, with values at edge centers.
cosϕᵈ = map(edges(s′)) do e
  (cos(point(s′, src(s′, e))[1]) + cos(point(s′, tgt(s′, e))[1])) / 2
end

α₀ = 0.354
α₂ = 0.25
α = map(point(s′)) do ϕ
  α₀ + α₂*((1/2)*(3*ϕ[1]^2 - 1))
end
A = 210
B = 2
f = 0.70
ρ = 1025
cw = 4186
H = 70
C = map(point(s′)) do ϕ
  f * ρ * cw * H
end
D = 0.6

#Tₛ₀ = map(point(s′)) do ϕ
#    12 .- 40*((1/2)*(3*(sin(ϕ[1]))^2 - 1))
#end
# Isothermal initial conditions:
Tₛ₀ = map(point(s′)) do ϕ
  15
end

# Store these values to be passed to the solver.

u₀ = ComponentArray{Float64}(Tₛ = Tₛ₀)

constants_and_parameters = (
  absorbed_radiation_α = α,
  outgoing_radiation_A = A,
  outgoing_radiation_B = B,
  energy_C = C,
  diffusion_D = D,
  cosϕᵖ = cosϕᵖ,
  diffusion_cosϕᵈ = cosϕᵈ)

#############################################
# Define how symbols map to Julia functions #
#############################################

# In this example, all operators come from the Discrete Exterior Calculus module
# from CombinatorialSpaces.
function generate(sd, my_symbol; hodge=GeometricHodge()) end

#######################
# Generate simulation #
#######################

sim = eval(gensim(budyko_sellers, dimension=1))
fₘ = sim(s, generate)

##################
# Run simulation #
##################

tₑ = 1e9

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")
@save "budyko_sellers.jld2" soln

soln = solve(prob, FBDF())

#############
# Visualize #
#############

lines(map(x -> x[1], point(s′)), soln(0.0).Tₛ)
lines(map(x -> x[1], point(s′)), soln(tₑ).Tₛ)

# Initial frame
frames = 100
fig = Figure(resolution = (800, 800))
ax1 = CairoMakie.Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(soln(tₑ).Tₛ))
Label(fig[1,1,Top()], "Surface temperature, Tₛ, [C°]")
Label(fig[2,1,Top()], "Line plot of temperature from North to South pole, every $(tₑ/frames) time units")

# Animation
record(fig, "budyko_sellers.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), soln(t).Tₛ)
end
