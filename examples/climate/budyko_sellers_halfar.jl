# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes

# External Dependencies
using MLStyle
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using CairoMakie
using SparseArrays
using GeometryBasics: Point2
Point2D = Point2{Float64};

halfar_eq2 = @decapode begin
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ  * d(h) ∧ (mag(♯(d(h)))^(n-1)) ∧ (h^(n+2)))
end
glens_law = @decapode begin
  Γ::Form1
  A::Form1
  (ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end
ice_dynamics_composition_diagram = @relation () begin
  dynamics(Γ,n)
  stress(Γ,n)
end
ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:Γ,:n]),
   Open(glens_law, [:Γ,:n])])
halfar = apex(ice_dynamics_cospan)
to_graphviz(halfar, verbose=false)

energy_balance = @decapode begin
  (Tₛ, ASR, OLR, HT)::Form0
  (C)::Constant

  Tₛ̇ == ∂ₜ(Tₛ) 

  Tₛ̇ == (ASR - OLR + HT) ./ C
end

absorbed_shortwave_radiation = @decapode begin
  (Q, ASR)::Form0
  α::Constant

  ASR == (1 .- α) .* Q
end

outgoing_longwave_radiation = @decapode begin
  (Tₛ, OLR)::Form0
  (A,B)::Constant

  OLR == A .+ (B .* Tₛ)
end

heat_transfer = @decapode begin
  (HT, Tₛ)::Form0
  (D,cosϕᵖ,cosϕᵈ)::Constant

  HT == (D ./ cosϕᵖ) .* ⋆(d(cosϕᵈ .* ⋆(d(Tₛ))))
end

insolation = @decapode begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end

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
  [Open(energy_balance,               [:Tₛ, :ASR, :OLR, :HT]),
   Open(absorbed_shortwave_radiation, [:Q, :ASR]),
   Open(outgoing_longwave_radiation,  [:Tₛ, :OLR]),
   Open(heat_transfer,                [:Tₛ, :HT, :cosϕᵖ]),
   Open(insolation,                   [:Q, :cosϕᵖ])])

budyko_sellers = apex(budyko_sellers_cospan)
write_json_acset(budyko_sellers, "budyko_sellers.json") # Save this Decapode as a JSON file
to_graphviz(budyko_sellers, verbose=false)


warming = @decapode begin
  Tₛ::Form0
  A::Form1

  A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.65e7)

end
to_graphviz(warming)

budyko_sellers_halfar_composition_diagram = @relation () begin
  budyko_sellers(Tₛ)

  warming(A, Tₛ)

  halfar(A)
end
to_graphviz(budyko_sellers_halfar_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")

budyko_sellers_halfar_cospan = oapply(budyko_sellers_halfar_composition_diagram,
  [Open(budyko_sellers, [:Tₛ]),
   Open(warming,        [:A, :Tₛ]),
   Open(halfar,         [:stress_A])])
budyko_sellers_halfar = apex(budyko_sellers_halfar_cospan)
to_graphviz(budyko_sellers_halfar)

budyko_sellers_halfar = expand_operators(budyko_sellers_halfar)
infer_types!(budyko_sellers_halfar, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(budyko_sellers_halfar, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(budyko_sellers_halfar)

s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s′, 100, point=Point2D.(range(-π/2 + π/32, π/2 - π/32, length=100), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())

cosϕᵖ = map(x -> cos(x[1]), point(s′))
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

Tₛ₀ = map(point(s′)) do ϕ
  15
end

lines(map(x -> x[1], point(s′)), Tₛ₀)

n = 3
ρ = 910
g = 9.8

h₀ = map(point(s′)) do (x,_)
  (((x)^2)+2.5) / 1e3
end

lines(map(x -> x[1], point(s′)), h₀)

# Store these values to be passed to the solver.
u₀ = ComponentArray(Tₛ=Tₛ₀, halfar_dynamics_h=h₀)

constants_and_parameters = (
  budyko_sellers_absorbed_radiation_α = α,
  budyko_sellers_outgoing_radiation_A = A,
  budyko_sellers_outgoing_radiation_B = B,
  budyko_sellers_energy_C = C,
  budyko_sellers_diffusion_D = D,
  budyko_sellers_cosϕᵖ = cosϕᵖ,
  budyko_sellers_diffusion_cosϕᵈ = cosϕᵈ,
  halfar_n = n,
  halfar_stress_ρ = ρ,
  halfar_stress_g = g)

# Symbols to functions

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :♯ => x -> begin
      e_vecs = map(edges(sd)) do e
        point(sd, sd[e, :∂v0]) - point(sd, sd[e, :∂v1])
      end
      neighbors = map(vertices(sd)) do v
        union(incident(sd, v, :∂v0), incident(sd, v, :∂v1))
      end
      n_vecs = map(neighbors) do es
        [e_vecs[e] for e in es]
      end
      map(neighbors, n_vecs) do es, nvs
        sum([nv*norm(nv)*x[e] for (e,nv) in zip(es,nvs)]) / sum(norm.(nvs))
      end
    end
    :mag => x -> begin
      norm.(x)
    end
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(s)
          append!(J, [s[e,:∂v0],s[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      avg_mat = sparse(I,J,V)
      avg_mat * x
    end
    :^ => (x,y) -> x .^ y
    :* => (x,y) -> x .* y
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

## Simulation generation

sim = eval(gensim(budyko_sellers_halfar, dimension=1))
fₘ = sim(s, generate)

tₑ = 1e6

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-4), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

@save "budyko_sellers_halfar.jld2" soln

## Visualize

lines(map(x -> x[1], point(s′)), soln(tₑ).Tₛ)

lines(map(x -> x[1], point(s′)), soln(tₑ).halfar_dynamics_h)

begin
frames = 100
fig = Figure(; size = (800, 800))
ax1 = CairoMakie.Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(soln(tₑ).Tₛ))
Label(fig[1,1,Top()], "Surface temperature, Tₛ, [C°]")
Label(fig[2,1,Top()], "Line plot of temperature from North to South pole, every $(tₑ/frames) time units")

record(fig, "budyko_sellers_halfar_T.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), soln(t).Tₛ)
end
end

begin
frames = 100
fig = Figure(; size = (800, 800))
ax1 = CairoMakie.Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(soln(tₑ).halfar_dynamics_h))
Label(fig[1,1,Top()], "Ice height, h")
Label(fig[2,1,Top()], "Line plot of ice height from North to South pole, every $(tₑ/frames) time units")

record(fig, "budyko_sellers_halfar_h.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), soln(t).halfar_dynamics_h)
end
end


