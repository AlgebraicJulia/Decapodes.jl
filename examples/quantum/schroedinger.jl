using Catlab, CombinatorialSpaces, Decapodes
using GLMakie, JLD2, LinearAlgebra, MultiScaleArrays, MLStyle, OrdinaryDiffEq
using GeometryBasics: Point2
Point2D = Point2{Float64}

Schroedinger = @decapode begin
  (i,h,m)::Constant
  V::Parameter
  Ψ::Form0
  
  ∂ₜ(Ψ) == ((-1 * (h^2)/(2*m))*Δ(Ψ) + V * Ψ) / (i*h)
end

infer_types!(Schroedinger, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(Schroedinger, op1_res_rules_1D, op2_res_rules_1D)
to_graphviz(Schroedinger)

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ => x -> begin
      lap = Δ(0,sd)
      lap * x
    end
    :^ => (x,y) -> x .^ y
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s′, 100, point=Point2D.(range(-1, 1, length=100), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())

sim = eval(gensim(Schroedinger, dimension=1))
fₘ = sim(s, generate)

Ψ = zeros(ComplexF64, nv(s))
Ψ[49] = 1e-16

u₀ = construct(PhysicsState, [VectorForm(Ψ)], ComplexF64[], [:Ψ])
constants_and_parameters = (
  i = im, # TODO: Relax assumption of Float64
  V = t -> begin
    z = zeros(ComplexF64, nv(s))
    z
  end,
  h = 6.5e-16, # Planck constant in [eV s]
  m = 5.49e-4, # mass of electron in [eV]
)

tₑ = 1e12
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

@save "schroedinger.jld2" soln

lines(map(x -> x[1], point(s′)), map(x -> x.re, findnode(soln(0.0), :Ψ)))
lines(map(x -> x[1], point(s′)), map(x -> x.re, findnode(soln(tₑ), :Ψ)))

begin
# Initial frame
frames = 100
fig = Figure(resolution = (800, 800))
ax1 = Axis(fig[1,1])
xlims!(ax1, extrema(map(x -> x[1], point(s′))))
ylims!(ax1, extrema(map(x -> x.re, findnode(soln(tₑ), :Ψ))))
Label(fig[1,1,Top()], "Ψ from Schroedinger Wave Equation")
Label(fig[2,1,Top()], "Line plot of real portion of Ψ, every $(tₑ/frames) time units")

# Animation
record(fig, "schroedinger.gif", range(0.0, tₑ; length=frames); framerate = 15) do t
  lines!(fig[1,1], map(x -> x[1], point(s′)), map(x -> x.re, findnode(soln(t), :Ψ)))
end
end
