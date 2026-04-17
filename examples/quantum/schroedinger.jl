# TODO: This simulation of the Schroedinger
# equation should be vetted by an outside reviewer.
using Catlab, CombinatorialSpaces, Decapodes, DiagrammaticEquations
using GLMakie, JLD2, LinearAlgebra, OrdinaryDiffEqTsit5
using ComponentArrays

Schroedinger = @decapode begin
  (i,h,m)::Constant
  V::Parameter
  Ψ::Form0
  
  ∂ₜ(Ψ) == ((-1 * (h^2)/(2*m))*Δ(Ψ) + V * Ψ) / (i*h)
end

function circle(n, c, float_type)
  s = EmbeddedDeltaSet1D{Bool, Point2d}()
  map(range(0, 2pi - (pi/(2^(n-1))); step=pi/(2^(n-1)))) do t
    add_vertex!(s, point=Point2d(cos(t),sin(t))*(c/2pi))
  end
  add_edges!(s, 1:(nv(s)-1), 2:nv(s))
  add_edge!(s, nv(s), 1)
  sd = EmbeddedDeltaDualComplex1D{Bool, float_type, Point2d}(s)
  subdivide_duals!(sd, Circumcenter())
  s, sd
end
s′, s = circle(7, 1, Float64);

sim = evalsim(Schroedinger, dimension=1, stateeltype=ComplexF64, preallocate=false)
fₘ = sim(s, nothing)

# Map the circle of unit circumference to the unit interval for setting ICs.
x_coords = [0, accumulate(+, s[:length])[1:end-1]...];
x₀, σ, k = 0.5, 0.08, 35.0
Ψ = map(x_coords) do x
  exp(-((x - x₀)^2) / (2σ^2)) * cis(k * x)
end;
∫squared_modulus(x) = sum(abs2.(x))
@show ∫squared_modulus(Ψ)
Ψ ./= √(∫squared_modulus(Ψ));
@show ∫squared_modulus(Ψ)

u₀ = ComponentArray(Ψ=Ψ);
constants_and_parameters = (
  i = im,
  V = t -> zeros(ComplexF64, nv(s)),
  h = 6.5e-16, # Planck constant in [eV s]
  m = 5.49e-4, # mass of electron in [eV]
);

tₑ = 1e10
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters);
soln = solve(prob, Tsit5(), dtmax=1e6);
@show ∫squared_modulus(soln(0.0).Ψ)
@show ∫squared_modulus(soln(tₑ).Ψ)
@show maximum(abs2.(soln(0.0).Ψ))
@show maximum(abs2.(soln(tₑ).Ψ))

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
ax1 = Makie.Axis(fig[1,1])
xlims!(ax1, extrema(x_coords))
ylims!(ax1, minimum(first, y_extrema), maximum(last, y_extrema))
Label(fig[1,1,Top()], "Ψ from Schroedinger Wave Equation")
Label(fig[2,1,Top()], "Probability amplitude |Ψ|², every $(tₑ/frames) time units")

# Animation
record(fig, "schroedinger.mp4", time_samples; framerate = 15) do t
  lines!(fig[1,1], x_coords, abs2.(soln(t).Ψ))
end
end
