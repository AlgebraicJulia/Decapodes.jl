# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes

# External Dependencies
using MLStyle
using MultiScaleArrays
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using GLMakie
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
budyko_sellers = @decapode begin
    (Q,Tₛ)::Form0
    (α,A,B,C,D,cosϕᵖ,cosϕᵈ)::Constant

    Tₛ̇ == ∂ₜ(Tₛ)
    ASR == (1 .- α) .* Q
    OLR == A .+ (B .* Tₛ)
    HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))

    Tₛ̇ == (ASR - OLR + HT) ./ C
end

# Infer the forms of dependent variables, and resolve which versions of DEC
# operators to use.
infer_types!(budyko_sellers, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(budyko_sellers, op1_res_rules_1D, op2_res_rules_1D)

to_graphviz(budyko_sellers)

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
Q = map(point(s′)) do ϕ
    450*cos(ϕ[1])
end

#Tₛ₀ = map(point(s′)) do ϕ
#    12 .- 40*((1/2)*(3*(sin(ϕ[1]))^2 - 1))
#end
# Isothermal:
Tₛ₀ = map(point(s′)) do ϕ
    15
end

# Store these values to be passed to the solver.
u₀ = construct(PhysicsState, [VectorForm(Q), VectorForm(Tₛ₀)], Float64[], [:Q, :Tₛ])
constants_and_parameters = (
    α = α,
    A = A,
    B = B,
    C = C,
    D = D,
    cosϕᵖ = cosϕᵖ,
    cosϕᵈ = cosϕᵈ)

#############################################
# Define how symbols map to Julia functions #
#############################################

hodge = GeometricHodge()
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :d₀ => x -> begin
      d₀ = d(s,0)
      d₀ * x
    end
    :dual_d₀ => x -> begin
      dual_d₀ = dual_derivative(s,0)
      dual_d₀ * x
    end
    :⋆₁ => x -> begin
      ⋆₁ = ⋆(s,1)
      ⋆₁ * x
    end
    :⋆₀⁻¹ => x -> begin
      ⋆₀⁻¹ = inv_hodge_star(s,0)
      ⋆₀⁻¹ * x
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

#######################
# Generate simulation #
#######################

sim = eval(gensim(budyko_sellers, dimension=1))
fₘ = sim(s, generate)

##################
# Run simulation #
##################

tₑ = 1e6

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

@save "budyko_sellers.jld2" soln

#############
# Visualize #
#############

lines(map(x -> x[1], point(s′)), findnode(soln(0.0), :Tₛ))
lines(map(x -> x[1], point(s′)), findnode(soln(tₑ), :Tₛ))

# Initial frame
fig = Figure(resolution = (800, 800))
ax1 = Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(findnode(soln(tₑ), :Tₛ)))
Label(fig[1,1,Top()], "Tₛ")

# Animation
frames = 100
record(fig, "budyko_sellers.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), findnode(soln(t), :Tₛ))
end
