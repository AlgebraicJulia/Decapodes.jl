using MultiScaleArrays
using OrdinaryDiffEq

struct VectorForm{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end

struct PhysicsState{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArray{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
    names::Vector{Symbol}
end

C = VectorForm(ones(Float64, 10))
V = VectorForm(ones(Float64, 100))

findname(u::PhysicsState, s::Symbol) = findfirst(isequal(s), u.names)
findnode(u::PhysicsState, s::Symbol) = u.nodes[findname(u, s)]



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

using JSON
using Distributions
using GLMakie

function plotform0(plot_mesh, c)
  fig, ax, ob = mesh(plot_mesh; color=c[point_map]);
  Colorbar(fig)
  ax.aspect = AxisAspect(3.0)
  fig
end

plot_mesh = parse_json_acset(EmbeddedDeltaSet2D{Bool, Point3{Float64}}, read("./docs/assets/meshes/plot_mesh.json", String))
periodic_mesh = parse_json_acset(EmbeddedDeltaDualComplex2D{Bool, Float64, Point3{Float64}}, read("./docs/assets/meshes/periodic_mesh.json", String));
point_map = JSON.parse(read("./docs/assets/meshes/point_map.json",String))

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
ddp = NamedDecapode(diffExpr)
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
plotform0(plot_mesh, findnode((soln(1)-u₀), :C))
plotform0(plot_mesh, findnode((soln(0.0000000000001)-u₀), :C))

times = range(0.0, tₑ, length=150)
colors = [findnode(soln(t), :C)[point_map] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "diff.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = findnode(soln(t), :C)[point_map]
end

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
advdiffdp = NamedDecapode(advdiff)
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
record(fig, "diff_adv.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = findnode(soln(t), :C)[point_map]
end

