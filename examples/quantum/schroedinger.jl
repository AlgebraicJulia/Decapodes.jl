using Catlab, CombinatorialSpaces, Decapodes
using GLMakie, JLD2, LinearAlgebra, OrdinaryDiffEq
using ComponentArrays
using GeometryBasics: Point2
Point2D = Point2{Float64}

Schroedinger = @decapode begin
  (i,h,m)::Constant
  V::Parameter
  Ψ::Form0
  
  ∂ₜ(Ψ) == ((-1 * (h^2)/(2*m))*Δ(Ψ) + V * Ψ) / (i*h)
end

s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s′, 100, point=Point2D.(range(-1, 1, length=100), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())

sim = eval(gensim(Schroedinger, dimension=1, stateeltype=ComplexF64))
fₘ = sim(s, nothing)

x_coords = map(x -> x[1], point(s′))
x₀, σ, k = -0.3, 0.08, 35.0
Ψ = map(x_coords) do x
  exp(-((x - x₀)^2) / (2σ^2)) * cis(k * x)
end

u₀ = ComponentArray(Ψ=Ψ)
constants_and_parameters = (
  i = im,
  V = t -> zeros(ComplexF64, nv(s)),
  h = 6.5e-16, # Planck constant in [eV s]
  m = 5.49e-4, # mass of electron in [eV]
)

tₑ = 1e12
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

@save "schroedinger.jld2" soln

lines(x_coords, abs2.(soln(0.0).Ψ))
lines(x_coords, abs2.(soln(tₑ).Ψ))

begin
# Initial frame
frames = 100
time_samples = range(0.0, tₑ; length=frames)
y_extrema = map(time_samples) do t
  extrema(abs2.(soln(t).Ψ))
end
fig = Figure(resolution = (800, 800))
ax1 = Axis(fig[1,1])
xlims!(ax1, extrema(x_coords))
ylims!(ax1, minimum(first, y_extrema), maximum(last, y_extrema))
Label(fig[1,1,Top()], "Ψ from Schroedinger Wave Equation")
Label(fig[2,1,Top()], "Probability amplitude |Ψ|², every $(tₑ/frames) time units")

# Animation
record(fig, "schroedinger.gif", time_samples; framerate = 15) do t
  lines!(fig[1,1], x_coords, abs2.(soln(t).Ψ))
end
end
