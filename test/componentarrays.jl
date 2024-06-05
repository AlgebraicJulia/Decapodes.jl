using ComponentArrays
using OrdinaryDiffEq
using GeometryBasics
using JSON
using Distributions

using Catlab
using Catlab.CategoricalAlgebra
using CombinatorialSpaces

using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using Test
using MLStyle
using LinearAlgebra

C = ones(Float64, 10)
V = ones(Float64, 100)

u₀ = ComponentArray(C=C,V=V)

dynamics(du, u, p, t) = begin
    du.C .= 0.1 * u.C
    return du
end
prob = ODEProblem(dynamics,u₀,(0,1))
soln = solve(prob, Tsit5())

function generate(sd, my_symbol)
  op = @match my_symbol begin
    :k => x->x/20
    _ => default_dec_generate(sd, my_symbol)
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
    (C, Ċ)::Form0{X}
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
fₘ = f(periodic_mesh, generate)
c_dist = MvNormal([5, 5], LinearAlgebra.Diagonal(map(abs2, [1.5, 1.5])))
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

u₀ = ComponentArray(C=c)
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

function closest_point(p1, p2, dims)
  p_res = collect(p2)
  for i in 1:length(dims)
    if dims[i] != Inf
      p = p1[i] - p2[i]
      f, n = modf(p / dims[i])
      p_res[i] += dims[i] * n
      if abs(f) > 0.5
        p_res[i] += sign(f) * dims[i]
      end
    end
  end
  Point3{Float64}(p_res...)
end

function flat_op(s::AbstractDeltaDualComplex2D, X::AbstractVector; dims=[Inf, Inf, Inf])
  # XXX: Creating this lookup table shouldn't be necessary. Of course, we could
  # index `tri_center` but that shouldn't be necessary either. Rather, we should
  # loop over incident triangles instead of the elementary duals, which just
  # happens to be inconvenient.
  tri_map = Dict{Int,Int}(triangle_center(s,t) => t for t in triangles(s))

  map(edges(s)) do e
  p = closest_point(point(s, tgt(s,e)), point(s, src(s,e)), dims)
  e_vec = (point(s, tgt(s,e)) - p) * sign(1,s,e)
  dual_edges = elementary_duals(1,s,e)
  dual_lengths = dual_volume(1, s, dual_edges)
  mapreduce(+, dual_edges, dual_lengths) do dual_e, dual_length
    X_vec = X[tri_map[s[dual_e, :D_∂v0]]]
    dual_length * dot(X_vec, e_vec)
  end / sum(dual_lengths)
  end
end

AdvDiff = quote
    (C, Ċ)::Form0{X}
    (V, ϕ, ϕ₁, ϕ₂)::Form1{X}

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
# Decapodes.compile(advdiffdp, [:C, :V])
# Decapodes.compile(expand_operators(advdiffdp), [:C, :V])
# gensim(expand_operators(advdiffdp), [:C, :V])

sim = eval(gensim(expand_operators(advdiffdp), [:C, :V]))
fₘ = sim(periodic_mesh, generate)
velocity(p) = [-0.5, -0.5, 0.0]
v = flat_op(periodic_mesh, DualVectorField(velocity.(periodic_mesh[triangle_center(periodic_mesh),:dual_point])); dims=[30, 10, Inf])
c_dist = MvNormal([7, 5], LinearAlgebra.Diagonal(map(abs2, [1.5, 1.5])))
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

u₀ = ComponentArray(C=c,V=v)
tₑ = 24
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

@test norm(soln.u[end].C - soln.u[1].C) >= 1e-4
@test norm(soln.u[end].V - soln.u[1].V) <= 1e-8

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

