""" diff_adv.jl

This file contains the code for the example presented in the Overview section
of the Decapodes.jl documentation. Comments have been removed, but description
of code is present in the Overview section.
"""

using Decapodes, Decapodes.Diagrams
using Catlab.Present, Catlab.Graphics

Variable = @decapode Decapodes2D begin
  C::Form0{X}
end;

draw_equation(decapode) = to_graphviz(decapode, node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"4.0"))

draw_equation(Variable)

TwoVariables = @decapode Decapodes2D begin
  C::Form0{X}
  dC::Form1{X}
end;

draw_equation(TwoVariables)

Equation = @decapode Decapodes2D begin
  C::Form0{X}
  dC::Form1{X}

  dC == d₀{X}(C)
end;

draw_equation(Equation)

@present DiffusionQuantities <: Decapodes2D begin
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
end;

Diffusion = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
  # Diffusion equation
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end;

draw_equation(Diffusion)

using Catlab.CategoricalAlgebra
using CombinatorialSpaces, CombinatorialSpaces.DiscreteExteriorCalculus
using CairoMakie

plot_mesh = loadmesh(Rectangle_30x10())
periodic_mesh = loadmesh(Torus_30x10())
point_map = loadmesh(Point_Map())

fig, ax, ob = wireframe(plot_mesh)
ax.aspect = AxisAspect(3.0)
fig

using Decapodes.Schedules

explicit_ts = diag2dwd(Diffusion)
to_graphviz(explicit_ts, orientation=LeftToRight)

using LinearAlgebra
using Decapodes.Examples, Decapodes.Simulations

funcs = sym2func(periodic_mesh)

funcs[:k] = Dict(:operator => 0.05 * I(ne(periodic_mesh)), :type => MatrixFunc())
funcs[:⋆₁] = Dict(:operator => ⋆(Val{1}, periodic_mesh, hodge=DiagonalHodge()),
                  :type => MatrixFunc());

func, code = gen_sim(explicit_ts, funcs, periodic_mesh; autodiff=false);

using Distributions
c_dist = MvNormal([7, 5], [1.5, 1.5])
c = [pdf(c_dist, [p[1], p[2]]) for p in periodic_mesh[:point]]

fig, ax, ob = mesh(plot_mesh; color=c[point_map])
ax.aspect = AxisAspect(3.0)
fig

using OrdinaryDiffEq

prob = ODEProblem(func, c, (0.0, 100.0))
sol = solve(prob, Tsit5());

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t)[point_map] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "diffusion.gif", range(0.0, 100.0; length=150); framerate = 30) do t
ob.color = sol(t)[point_map]
end

Diffusion = @decapode DiffusionQuantities begin
  C::Form0{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(C))
end

Advection = @decapode DiffusionQuantities begin
  C::Form0{X}
  (V, ϕ)::Form1{X}
  ϕ == ∧₀₁{X}(C,V)
end

Superposition = @decapode DiffusionQuantities begin
  (C, Ċ)::Form0{X}
  (ϕ, ϕ₁, ϕ₂)::Form1{X}

  ϕ == ϕ₁ + ϕ₂
  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(C) == Ċ
end

draw_equation(Diffusion)

draw_equation(Advection)

draw_equation(Superposition)
using Catlab.Programs

compose_diff_adv = @relation (C, V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end

to_graphviz(compose_diff_adv, box_labels=:name, junction_labels=:variable,
            graph_attrs=Dict(:start => "2"))

using Decapodes.OpenDiagrams
DiffusionAdvection = oapply(compose_diff_adv,
                  [OpenDiagram(Diffusion, [:C, :ϕ]),
                   OpenDiagram(Advection, [:C, :ϕ, :V]),
                   OpenDiagram(Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])])

draw_equation(DiffusionAdvection.functor)
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

explicit_ts = diag2dwd(DiffusionAdvection.functor)
to_graphviz(explicit_ts, orientation=LeftToRight)

using CombinatorialSpaces.DiscreteExteriorCalculus: ∧
funcs[:∧₀₁] = Dict(:operator => (r, c,v)->r .= ∧(Tuple{0,1}, periodic_mesh, c, v), :type => InPlaceFunc())

func, code = gen_sim(explicit_ts, funcs, periodic_mesh; autodiff=false, params = [:V]);

velocity(p) = [-0.5, -0.5, 0.0]
v = flat_op(periodic_mesh, DualVectorField(velocity.(periodic_mesh[triangle_center(periodic_mesh),:dual_point])); dims=[30, 10, Inf])

prob = ODEProblem(func, c, (0.0, 100.0))
sol = solve(prob, Tsit5(), p=v);

# Plot the result
times = range(0.0, 100.0, length=150)
colors = [sol(t)[point_map] for t in times]

# Initial frame
fig, ax, ob = mesh(plot_mesh, color=colors[1], colorrange = extrema(vcat(colors...)))
ax.aspect = AxisAspect(3.0)
Colorbar(fig[1,2], ob)
framerate = 30

# Animation
record(fig, "diff_adv.gif", range(0.0, 100.0; length=150); framerate = 30) do t
ob.color = sol(t)[point_map]
end
