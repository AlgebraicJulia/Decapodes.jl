# We use here the formulation studied by Jordan, Kinderlehrer, and Otto in "The
# Variational Formulation of the Fokker-Planck Equation" (1996).

# The formualation they studied is that where the drift coefficient is the
# gradient of (some potential) Ψ.

# Load libraries.
using Catlab, CombinatorialSpaces, Decapodes
using GLMakie, LinearAlgebra, MLStyle, MultiScaleArrays, OrdinaryDiffEq

# Specify physics.
Fokker_Planck = @decapode begin
  (ρ,Ψ)::Form0
  β⁻¹::Constant
  ∂ₜ(ρ) == ∘(⋆,d,⋆)(d(Ψ)∧ρ) + β⁻¹*Δ(ρ)
end

Fokker_Planck = expand_operators(Fokker_Planck)
infer_types!(Fokker_Planck)
resolve_overloads!(Fokker_Planck)

# Specify domain.
include("../grid_meshes.jl")
include("examples/grid_meshes.jl")
s′ = triangulated_grid(1,1,0.05,0.05,Point3D)
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())

# Specify initial conditions.
Ψ = map(p -> √(2)/2 - √((p[1]-.5)^2 + (p[2]-.5)^2) , point(s))
ρ = fill(1 / sum(s[:area]), nv(s))
constants_and_parameters = (β⁻¹ = 0.001,)
u₀ = construct(PhysicsState, [VectorForm(Ψ), VectorForm(ρ)], Float64[], [:Ψ, :ρ])

# Compile.
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    _ => default_dec_generate(sd, my_symbol)
  end
  return (args...) ->  op(args...)
end
sim = eval(gensim(Fokker_Planck))
fₘ = sim(s, generate)

# Run.
tₑ = 1.0
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())

# Save solution data.
@save "fokker_planck.jld2" soln

# Create GIF
begin
frames = 800
# Initial frame
fig = CairoMakie.Figure(resolution = (800, 800))
p1 = CairoMakie.mesh(fig[1,1], s′, color=findnode(soln(0), :ρ), colormap=:jet, colorrange=extrema(findnode(soln(0), :ρ)))
Colorbar(fig[1,2], colormap=:jet, colorrange=extrema(findnode(soln(0), :U)))

# Animation
record(fig, "fokker_planck.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    p1.plot.color = findnode(soln(t), :ρ)
end

end

