using MultiScaleArrays
using OrdinaryDiffEq
using GeometryBasics
using JSON
using Distributions
using Artifacts
# using GLMakie

using Catlab
using Catlab.CategoricalAlgebra
using CombinatorialSpaces
using Decapodes
using Test
using MLStyle
using LinearAlgebra

C = VectorForm(ones(Float64, 10))
V = VectorForm(ones(Float64, 100))

u₀ = construct(PhysicsState, [C,V], Float64[], [:C, :V])
@test length(findnode(u₀, :C)) == 10
@test length(findnode(u₀, :V)) == 100

dynamics(du, u, p, t) = begin
    findnode(du, :C).values .= 0.1 * findnode(u, :C).values
    return du
end
prob = ODEProblem(dynamics,u₀,(0,1))
soln = solve(prob, Tsit5())

function generate(sd, my_symbol)
  op = @match my_symbol begin
    :k => x->x/20
    :⋆₀ => x->⋆(0,sd,hodge=DiagonalHodge())*x
    :⋆₁ => x->⋆(1, sd, hodge=DiagonalHodge())*x
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=DiagonalHodge())
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=DiagonalHodge())*x
    :d₀ => x->d(0,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
    :plus => (+)
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length(args[1]))"); op(args...);end
  return (args...) ->  op(args...)
end


plot_mesh = loadmesh(Rectangle_30x10())
periodic_mesh = loadmesh(Torus_30x10())
point_map = loadmesh(Point_Map())

# function plotform0(plot_mesh, c)
#   fig, ax, ob = mesh(plot_mesh; color=c[point_map]);
#   Colorbar(fig)
#   ax.aspect = AxisAspect(3.0)
#   fig
# end

DiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    ϕ::Form1{X}

    # Fick's first law
    ϕ ==  ∘(d₀, k)(C)
    # Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

diffExpr = parse_decapode(DiffusionExprBody)
ddp = SummationDecapode(diffExpr)
gensim(expand_operators(ddp), [:C])
f = eval(gensim(expand_operators(ddp), [:C]))
fₘ = f(periodic_mesh)
c_dist = MvNormal([5, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

u₀ = construct(PhysicsState, [VectorForm(c)],Float64[], [:C])
tₑ = 10
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

soln(0.9)
# plotform0(plot_mesh, findnode((soln(1)-u₀), :C))
# plotform0(plot_mesh, findnode((soln(0.0000000000001)-u₀), :C))

# times = range(0.0, tₑ, length=150)
# colors = [findnode(soln(t), :C)[point_map] for t in times]

# Initial frame
# fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
# ax.aspect = AxisAspect(3.0)
# Colorbar(fig[1,2], ob)
# framerate = 30

# Animation
# record(fig, "diff.gif", range(0.0, tₑ; length=150); framerate = 30) do t
#     ob.color = findnode(soln(t), :C)[point_map]
# end

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

advdiff = parse_decapode(AdvDiff)
advdiffdp = SummationDecapode(advdiff)
Decapodes.compile(advdiffdp, [:C, :V])
Decapodes.compile(expand_operators(advdiffdp), [:C, :V])
gensim(expand_operators(advdiffdp), [:C, :V])
sim = eval(gensim(expand_operators(advdiffdp), [:C, :V]))
fₘ = sim(periodic_mesh)
velocity(p) = [-0.5, -0.5, 0.0]
v = flat_op(periodic_mesh, DualVectorField(velocity.(periodic_mesh[triangle_center(periodic_mesh),:dual_point])); dims=[30, 10, Inf])
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(v)],Float64[], [:C, :V])
tₑ = 24
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

@test norm(findnode(soln.u[end], :C) - findnode(soln.u[1], :C)) >= 1e-4
@test norm(findnode(soln.u[end], :V) - findnode(soln.u[1], :V)) <= 1e-8

# Plot the result
# times = range(0.0, tₑ, length=150)
# colors = [findnode(soln(t), :C)[point_map] for t in times]

# Initial frame
# fig, ax, ob = mesh(plot_mesh, color=colors[end], colorrange = extrema(vcat(colors...)))
# ax.aspect = AxisAspect(3.0)
# Colorbar(fig[1,2], ob)
# framerate = 30

# Animation
# record(fig, "diff_adv.gif", range(0.0, tₑ; length=150); framerate = 30) do t
#     ob.color = findnode(soln(t), :C)[point_map]
# end
