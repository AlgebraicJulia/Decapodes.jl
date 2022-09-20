using Distributions
using GLMakie
using MultiScaleArrays
using OrdinaryDiffEq

include("./spherical_meshes.jl")

struct VectorForm{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end

struct PhysicsState{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
    names::Vector{Symbol}
end

findname(u::PhysicsState, s::Symbol) = findfirst(isequal(s), u.names)
findnode(u::PhysicsState, s::Symbol) = u.nodes[findname(u, s)]

const EARTH_RADIUS = 6371
const MIN_THERMOSPHERE_ALTITUDE = 91
const WITH_NORTH_POLE = true
tie_gcm_grid, npi, spi = makeSphere(WITH_NORTH_POLE ? 0 : 5, 180, 5, 0, 360, 5,
                          EARTH_RADIUS+MIN_THERMOSPHERE_ALTITUDE)

c_dist = MvNormal([EARTH_RADIUS+MIN_THERMOSPHERE_ALTITUDE, 0, 0],
                  [EARTH_RADIUS/6, EARTH_RADIUS/6, EARTH_RADIUS/6])
c = map(p -> pdf(c_dist, p), tie_gcm_grid[:point])
# TODO: Set value equal to distance from z-axis.
v = collect(1:2522)
#velocity(p) = [-0.5, -0.5, 0.0]
#v = flat_op(periodic_mesh, DualVectorField(velocity.(periodic_mesh[triangle_center(periodic_mesh),:dual_point])); dims=[30, 10, Inf])

AdvDiff = quote
    C::Form0{X}
    Ċ::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}

    # Fick's first law
    ϕ₁ ==  (d₀∘k)(C)
    ϕ₂ == ∧₀₁(C,V)
    ϕ == ϕ₁ + ϕ₂
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁,⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

fₘ = (advdiff |> parse_decapode |> NamedDecapode |> expand_operators |>
      x -> gensim(x, [:C, :V]) |> eval)(tie_gcm_grid)

u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(v)],Float64[], [:C, :V])
tₑ = 24
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

@test norm(findnode(soln.u[end], :V) - findnode(soln.u[1], :V)) <= 1e-8

# Plot the result
times = range(0.0, tₑ, length=150)
colors = [findnode(soln(t), :C)[point_map] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[end], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "space_weather.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = findnode(soln(t), :C)[point_map]
end

