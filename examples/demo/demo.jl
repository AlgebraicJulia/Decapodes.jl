using Semagrams, Decapodes, Catlab, DifferentialEquations, CairoMakie, CombinatorialSpaces, LinearAlgebra

using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules
using Decapodes.Debug
using Decapodes.OpenDiagrams

using Logging: global_logger
using TerminalLoggers: TerminalLogger

using Distributions

coords = [[-0.0,0,0],[-0.0,4,0],[2.75,2,0],[3.5 + 2/sqrt(3),2,0],
          [3.5,4,0], [3.5,0,0], [7,0,0],  [3.5 + 2/sqrt(3)+2/sqrt(3),4,0],  [7,2,0]]

primal_s = EmbeddedDeltaSet2D{Bool,Point3{Float64}}()
add_vertices!(primal_s, 9, point=[Point3{Float64}(c) for c in coords])
glue_triangle!(primal_s, 1, 2, 3)
glue_triangle!(primal_s, 1, 3, 6)
glue_triangle!(primal_s, 2, 3, 5)
glue_triangle!(primal_s, 3, 6, 4)
glue_triangle!(primal_s, 3, 4, 5)
glue_triangle!(primal_s, 4, 5, 8)
glue_triangle!(primal_s, 4, 9, 8)
glue_triangle!(primal_s, 4, 7, 9)
glue_triangle!(primal_s, 6, 4, 7)
s = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3{Float64}}(primal_s)

subdivide_duals!(s, Circumcenter())
///

fig, ax, ob = wireframe(s)
wireframe!(primal_s, linecolor=:black)
fig
///

height = 10
width = 30
ncols = 30

get_index(s, p::Point3{Float64}; σ=1e-3) = findall(ps -> norm(ps - p) < σ, s[:point])

function add_unique_point(s, p::Point3{Float64}; σ=1e-3)
    pnts = get_index(s, p, σ=σ)
    if isempty(pnts)
        add_vertex!(s, point=p)
    else
        length(pnts) == 1 || error("Multiple points found, reduce σ range")
        only(pnts)
    end
end

primal_s = EmbeddedDeltaSet2D{Bool,Point3{Float64}}()
col_h = height / ncols

bound_pts = Point3{Float64}.([[0,1,0], [0.875,1,0], [0.6875, 0.5,0], [0,0,0], [0.875,0,0]] .* col_h)
bound_tris = [[1,2,3],[1,3,4],[3,4,5]]

for c in 1:ncols
    off_y = col_h  * (c-1)
    right_edge = Point3{Float64}([width, off_y + col_h, 0])
    left_edge = Point3{Float64}([0, off_y, 0])
    # Left boundary
    vts = [add_unique_point(primal_s, left_edge + bp) for bp in bound_pts]
    tris = [glue_triangle!(primal_s, vts[bt]...) for bt in bound_tris]
    
    # Right boundary
    vts = [add_unique_point(primal_s, right_edge - bp) for bp in bound_pts]
    tris = [glue_triangle!(primal_s, vts[bt]...) for bt in bound_tris]

    # Filling
    # Calculate width to get as close to equilateral triangles
    eff_width = width - 2 * 0.875 * col_h
    ideal_w = col_h / sqrt(3)
    n_tris = round(Int64, eff_width / ideal_w)
    
    lu_w = eff_width / n_tris
    c_w = (width - 2 * 0.6875 * col_h) / (n_tris+1)
    lup, cp = Point3{Float64}.([[lu_w, 0, 0],[c_w, 0,0]])
    l_vts = [add_unique_point(primal_s, left_edge + bound_pts[5] + lup * (i-1)) for i in 1:(n_tris+1)]
    c_vts = [add_unique_point(primal_s, left_edge + bound_pts[3] + cp * (i-1)) for i in 1:(n_tris+2)]
    u_vts = [add_unique_point(primal_s, left_edge + bound_pts[2] + lup * (i-1)) for i in 1:(n_tris+1)]
    for i in 0:(n_tris*2)
        ind = i ÷ 2 + 1
        if i % 2 == 0
            glue_triangle!(primal_s, l_vts[ind], c_vts[ind], c_vts[ind+1])
            glue_triangle!(primal_s, u_vts[ind], c_vts[ind], c_vts[ind+1])
        else
            glue_triangle!(primal_s, l_vts[ind], l_vts[ind+1], c_vts[ind+1])
            glue_triangle!(primal_s, u_vts[ind], u_vts[ind+1], c_vts[ind+1])
        end
    end
end
orient!(primal_s)
s = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3{Float64}}(primal_s)

subdivide_duals!(s, Circumcenter())
fig, ax, ob = wireframe(s)
ax.aspect = AxisAspect(width/height)
wireframe!(primal_s, color=:black)
fig
///

@present Flow2DQuantities <: Decapodes2D begin
  X::Space
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  avg₀₁::Hom(Form0(X), Form1(X))
  prod₁::Hom(Form1(X)⊗Form1(X), Form1(X))
end;
///

Diffusion = @decapode Flow2DQuantities begin
  (T, Ṫ)::Form0{X}
  V::Form1{X}
  ϕ::Form1{X}

  # Fick's first law
  ϕ ==  k(d₀{X}(T)) + ∧₀₁{X}(T,V)
  # Diffusion equation
  Ṫ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(ϕ)))
  ∂ₜ{Form0{X}}(T) == Ṫ
  ∂ₜ{Form1{X}}(V) == prod₁(avg₀₁(T), V)
end;
to_graphviz(Diffusion)
///

(dwd = diag2dwd(Diffusion)) |> to_graphviz
///

using SparseArrays
funcs = sym2func(s)

function avg_mat(::Type{Val{(0,1)}},s)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      append!(J, [s[e,:∂v0],s[e,:∂v1]])
      append!(I, [e,e])
      append!(V, [0.5, 0.5])
  end
  sparse(I,J,V)
end

funcs[:k] = Dict(:operator => 0.05 * I(ne(s)), :type => MatrixFunc())
funcs[:∧₀₁] = Dict(:operator => (r, c,v)->r .= ∧(Tuple{0,1}, s, c, v), :type => InPlaceFunc())
funcs[:prod₁] = Dict(:operator => (x′, x,y)->x′ .= x .* y, :type => InPlaceFunc())
funcs[:avg₀₁] = Dict(:operator => avg_mat(Val{(0,1)}, s), :type => MatrixFunc())

Examples.contract_matrices!(dwd, funcs)
to_graphviz(dwd)
///

c_dist = MvNormal([7, 5], [1.5, 1.5])
c_trail =  MvNormal([3,5],[1,1])
c = [pdf(c_dist, [p[1], p[2]]) + pdf(c_trail, [p[1], p[2]])*0.25  for p in primal_s[:point]]

velocity(p) = [-1.0, 0.0, 0.0]
v = ♭(s, DualVectorField(velocity.(s[triangle_center(s),:dual_point]))).data;

u0 = vcat(c, v)

func, code = gen_sim(dwd, funcs, s; autodiff=false);
///

prob = ODEProblem(func, u0, (0.0, 20.0))
sol = solve(prob, Tsit5(), progress=true, progress_steps=10, p=v);
///

fig, ax, ob = mesh(primal_s, color=sol(17.0)[1:nv(s)])
ax.aspect = AxisAspect(width/height)
fig
///

c_range = 1:nv(s)

times = range(0.0, 17.0, length=150)
colors = [sol(t)[c_range] for t in times]

fig, ax, ob = mesh(primal_s, color=colors[1], colorrange = (0.0, maximum(vcat(colors...))))
ax.aspect = AxisAspect(width/height)
Colorbar(fig[1,2], ob)
framerate = 30

record(fig, "flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
ob.color = colors[i]
end
///
