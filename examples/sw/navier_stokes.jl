using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
using SparseArrays
using GLMakie
using Logging
# using CairoMakie 

using GeometryBasics: Point3
Point3D = Point3{Float64}

# Begin Navier Stokes

# Navier Stokes example
DiffusionExprBody = quote
  (T, Ṫ)::Form0{X}
  ϕ::DualForm1{X}
  k::Constant{X}
  # Fick's first law
  ϕ ==  ⋆₁(k*d₀(T))
  # Diffusion equation
  Ṫ == ⋆₀⁻¹(dual_d₁(ϕ))
end

Diffusion = SummationDecapode(parse_decapode(DiffusionExprBody))
to_graphviz(Diffusion)
to_graphviz(Diffusion, graph_attrs=Dict(:rankdir => "LR"))

AdvectionExprBody = quote
  #(M,V)::Form1{X}  #  M = ρV
  V::Form1{X}  #  M = ρV
  (ρ, P, T, Ṫ)::Form0{X}
  #V == M/avg₀₁(ρ)
  #ρ == P / R₀(T)
  ρ == P ./ R₀(T)
  Ṫ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(T))))
end

Advection = SummationDecapode(parse_decapode(AdvectionExprBody))
to_graphviz(Advection)
to_graphviz(Advection, graph_attrs=Dict(:rankdir => "LR"))

SuperpositionExprBody = quote
  (T, Ṫ, Ṫ₁, Ṫₐ)::Form0{X}
  Ṫ == Ṫ₁ + Ṫₐ
  ∂ₜ(T) == Ṫ 
end
Superposition = SummationDecapode(parse_decapode(SuperpositionExprBody))
to_graphviz(Superposition)
to_graphviz(Superposition, graph_attrs=Dict(:rankdir => "LR"))

compose_continuity = @relation () begin
  diffusion(T, Ṫ₁)
  #advection(M, ρ, T, Ṫₐ)
  advection(ρ, T, Ṫₐ)
  superposition(T, Ṫ, Ṫ₁, Ṫₐ)
end
to_graphviz(compose_continuity, junction_labels=:variable, box_labels=:name, prog="circo")

continuity_cospan = oapply(compose_continuity,
                [Open(Diffusion, [:T, :Ṫ]),
                 Open(Advection, [:ρ, :T, :Ṫ]),
                 Open(Superposition, [:T, :Ṫ, :Ṫ₁, :Ṫₐ])])

continuity = apex(continuity_cospan)
to_graphviz(continuity)
to_graphviz(continuity, graph_attrs=Dict(:rankdir => "LR"))

NavierStokesExprBody = quote
  #(M, Ṁ, G, V)::Form1{X}
  (M, Ṁ, V)::Form1{X}
  (ρ, P, Ṗ)::Form0{X}
  B::Form2{X}
  (negone,two,three,c,kᵥ,q,nₑ)::Constant{X}
  #V == M/avg₀₁(ρ)
  #M == V*avg₀₁(ρ)
  #M == V .* avg₀₁(ρ)
  #Ṁ == neg₁(L₁′(V, V))*avg₀₁(ρ) + 
  Ṁ == neg₁(L₁′(V, V)) .* avg₀₁(ρ) + 
        #kᵥ*(Δ₁(V) + d₀(δ₁(V))/three) +
        kᵥ*(Δ₁(V) + d₀(δ₁(V))/three) .* avg₀₁(ρ) +
        #d₀(i₁′(V, V)/two)*avg₀₁(ρ) +
        d₀(i₁′(V, V)/two) .* avg₀₁(ρ) +
        neg₁(d₀(P)) #+ # Possibly wrong
        #G*avg₀₁(ρ) + 
        # (q*V cross product B) # 1 form in primal mesh with 0 form in dual mesh
        #cp_2_1(B, q*⋆₁⁻¹(V)) * negone
        #cp_2_1(B, (nₑ*q)*⋆₁(V)) * negone
        #∧₁₀′(B, (nₑ*q)*⋆₁(V)) * negone
        #∧₁₀′((nₑ*q)*⋆₁(V), B) * negone
        #∧₁₀′((nₑ*q/c)*⋆₁(V), B) * negone
  ∂ₜ(M) == Ṁ
  Ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(P)))) # *Lie(3Form) = Div(*3Form x v) --> conservation of pressure
  ∂ₜ(P) == Ṗ
end

NavierStokes = SummationDecapode(parse_decapode(NavierStokesExprBody))
to_graphviz(NavierStokes)
to_graphviz(NavierStokes, graph_attrs=Dict(:rankdir => "LR"))

compose_heatXfer = @relation () begin
  continuity(V, ρ, P)
  navierstokes(V, ρ, P)
end
to_graphviz(compose_heatXfer, junction_labels=:variable, box_labels=:name, prog="circo")

heatXfer_cospan = oapply(compose_heatXfer,
                [Open(continuity, [:advection_V, :ρ, :advection_P]),
                 Open(NavierStokes, [:V, :ρ, :P])])

HeatXfer = apex(heatXfer_cospan)
to_graphviz(HeatXfer)
to_graphviz(HeatXfer, graph_attrs=Dict(:rankdir => "LR"))

# x is a Form1, y is a DualForm0 or Form2
# Returns a Form1.
function cp_2_1(x,y)
  #x_cache = t2c * x
  #y_cache = e2c * y
  x_cache = e2c * x
  y_cache = t2c * y
  broadcast!(*, y_cache, x_cache, y_cache)
  cross * y_cache
end

function edge_to_support(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      de = elementary_duals(Val{1},s, e)
      ev = CombinatorialSpaces.volume(Val{1}, s, e)
      dv = sum([dual_volume(Val{1}, s, d) for d in de])
      for d in de
          dt = incident(s, d, :D_∂e0)
          append!(I, dt)
          append!(J, fill(e, length(dt)))
          append!(V, fill(1/(dv*ev), length(dt)))
      end
  end
  sparse(I,J,V)
end

function support_to_tri(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:nv(s)
      dt = elementary_duals(Val{0},s, t)
      for d in dt
          push!(I, t)
          push!(J, d)
          push!(V, 1)
      end
  end
  sparse(I,J,V)
end

diag_vols(s) = spdiagm([dual_volume(Val{2}, s, dt) for dt in 1:nparts(s, :DualTri)])

wedge_mat(::Type{Val{(1,1)}}, s) = 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))

function pd_wedge(::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), kw...)
  wedge_t[(1,1)] * broadcast(*, α, β)
end

#∧₁₀′(⋆₁(v,earth,hodge), Br_temp)
#∧₁₀′(Br_temp, ⋆₁(v,earth,hodge))
#
#∧₁₀′(v, ⋆₂(d₁(v,earth),earth,hodge))
#⋆₂(d₁(v,earth),earth,hodge)
#∧₁₀′(⋆₂(d₁(v,earth),earth,hodge),v)
#i₀′(v,d₁(v,earth),earth,hodge)
##∧₁₀′(⋆₂(d₁(v,earth),earth,hodge), v)
#∧₁₀′(v, ⋆₂(d₁(v,earth),earth,hodge))
#i₀′(v,d₁(v,earth),earth,hodge) + d₀(i₁′(v,v,earth,hodge),earth)

hodge = GeometricHodge()
d₀(x,sd) = d(0,sd)*x
d₁(x,sd) = d(1,sd)*x
⋆₀(x,sd,hodge) = ⋆(0,sd,hodge=hodge)*x
⋆₁(x,sd,hodge) = ⋆(1,sd,hodge=hodge)*x
⋆₂(x,sd,hodge) = ⋆(2,sd,hodge=hodge)*x
⋆₀⁻¹(x,sd,hodge) = inv_hodge_star(0,sd,hodge=hodge)*x
⋆₁⁻¹(x,sd,hodge) = inv_hodge_star(1,sd,hodge=hodge)*x
⋆₂⁻¹(x,sd,hodge) = inv_hodge_star(2,sd,hodge=hodge)*x
∧₁₀′(x,y) = cp_2_1(x,y)
∧₁₁′(x,y,sd) = pd_wedge(Val{(1,1)},sd, x, y)
i₀′(x,y,sd,hodge) = -1.0 * (∧₁₀′(x, ⋆₂(y,sd,hodge)))
#i₀′(x,y,sd,hodge) = -1.0 * (∧₁₀′(⋆₂(y,sd,hodge), x))
i₁′(x,y,sd,hodge) = ⋆₀⁻¹(∧₁₁′(x, ⋆₁(y,sd,hodge), sd),sd,hodge)
L₁′(x,y,sd,hodge) = i₀′(x,d₁(y,sd),sd,hodge) + d₀(i₁′(x,y,sd,hodge),sd)
boltzmann_constant = 1.38064852e-23
mol_mass = 28.96
density = 0.000210322
kₜ = 0.0246295028571 # Thermal conductivity
cₚ = 1004.703 # Specific Heat at constant pressure
k₁ = kₜ / (density * cₚ) # Heat diffusion constant in fluid

#vect(s, e) = (s[s[e,:∂v1], :point] - s[s[e,:∂v0], :point]) * sign(1, s, e)
#vect(s, e::AbstractVector) = [vect(s, el) for el in e]
#t_vects(s,t) = vect(s, triangle_edges(s,t)) .* ([1,-1,1] * sign(2,s,t))
# These helper functions convert an edge to a vector, multiplying by the sign of
# a simplex where appropriate.
edge_to_vector(s, e) = (s[e, [:∂v1, :point]] - s[e, [:∂v0, :point]]) * sign(1, s, e)
edge_to_vector(s, e::AbstractVector) = [edge_to_vector(s, el) for el in e]
t_vects(s,t) = edge_to_vector(s, triangle_edges(s,t)) .* ([1,-1,1] * sign(2,s,t))
# Return a dictionary where a key is a (triangle_idx, edge_idx) a pair, and a
# value is an index for that pair (from 1 to num_triangles*3).
# The purpose of this dictionary is so we can store values on an edge with
# respect to the triangle(s) it is a part of individually.
function comp_support(sd)
  vects = []
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(vects, (t, inds[i]))
      end
  end
  v2comp = Dict{Tuple{Int64, Int64}, Int64}()
  for (i, v) in enumerate(vects)
      v2comp[v] = i
  end
  v2comp
end
v2comp = comp_support(earth)

# Return a (num_triangles*3) by (num_edges) sparse matrix.
# Row index is the index of the (triangle_idx, edge_idx) pair according to
# comp_support.
# Column index is the index of an edge.
# Values are 1/length of the corresponding edge.
# Multiplying by this matrix thus normalizes by edge length.
function edge2comp(sd, v2comp)
  I = Vector{Int64}(); J = Vector{Int64}(); V = Vector{Float64}()
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(I, v2comp[(t,inds[i])])
          push!(J, inds[i])
          #push!(V, 1 / volume(Val{1}, sd, inds[i]))
          push!(V, 1 / SimplicialSets.volume(Val{1}, sd, inds[i]))
      end
  end
  sparse(I,J,V)
end
e2c = edge2comp(earth, v2comp)

# Return a (num_triangles*3) by (num_triangles) sparse matrix.
# Row index is the index of the (triangle_idx, edge_idx) pair according to
# comp_support.
# Column index is the index of a triangle.
# Values are 1.
function tri2comp(sd, v2comp)
  I = Vector{Int64}(); J = Vector{Int64}(); V = Vector{Float64}()
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(I, v2comp[(t,inds[i])])
          push!(J, t)
          push!(V, 1)
      end
  end
  sparse(I,J,V)
end
t2c = tri2comp(earth, v2comp)

# Return a (num_edges) by (num_triangles*3) sparse matrix.
function changes(s, v2comp)
  orient_vals = [1,-1,1]
  I = Vector{Int64}(); J = Vector{Int64}(); V = Vector{Float64}()
  for t in 1:ntriangles(s)
    inds = triangle_edges(s, t)
    e_vects = t_vects(s,t)
    for i in 1:3
      ns = [(i+1)%3 + 1, i%3+1]
      ort = e_vects[i] × (e_vects[i] × e_vects[ns[1]])
      n_ort = normalize(ort)
      append!(I, inds[ns[1]])
      append!(J, v2comp[(t,inds[i])])
      append!(V, dot(n_ort, e_vects[ns[1]]) * orient_vals[ns[1]] * sign(1, s, ns[1])* orient_vals[i]* sign(2,s,t) / 3.0)
      append!(I, inds[ns[2]])
      append!(J, v2comp[(t,inds[i])])
      append!(V, dot(n_ort, e_vects[ns[2]]) * orient_vals[ns[2]] * sign(1, s, ns[2])* orient_vals[i]* sign(2,s,t) / 3.0)
    end
  end
  sparse(I,J,V, ne(s), ntriangles(s)*3)
end
cross = changes(earth, v2comp)

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

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    #:k => x->2000x
    :μ => x->-0.0001x
    # :μ => x->-2000x
    :β => x->2000*x
    :γ => x->1*x
    :⋆₀ => x->⋆(0,sd,hodge=hodge)*x
    :⋆₁ => x->⋆(1, sd, hodge=hodge)*x
    #:⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=hodge)
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, hodge=hodge)*x
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=hodge)*x
    :d₀ => x->d(0,sd)*x
    :d₁ => x->d(1,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
    :∧₁₀′ => (x,y) -> cp_2_1(x,y)
    #:cp_2_1 => (x,y)-> begin
    #  x_cache = t2c * x
    #  y_cache = e2c * y
    #  broadcast!(*, y_cache, x_cache, y_cache)
    #  cross * y_cache
    #end
    :plus => (+)
    :(-) => x-> -x
    # :L₀ => (v,x)->dual_derivative(1, sd)*(i0(v, x))
    :L₀ => (v,x)->begin
      # dual_derivative(1, sd)*⋆(1, sd)*wedge_product(Tuple{1,0}, sd, v, x)
      dual_derivative(1, sd)*⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x)
    end
    :i₀ => i0 
    :Δ₀ => x -> begin # dδ
      δ(1, sd, d(0, sd)*x, hodge=hodge) end
    # :Δ₀ => x -> begin # d ⋆ d̃ ⋆⁻¹
    #   y = dual_derivative(1,sd)*⋆(1, sd, hodge=hodge)*d(0,sd)*x
    #   inv_hodge_star(0,sd, y; hodge=hodge)
    # end
    :Δ₁ => x -> begin # dδ + δd
      δ(2, sd, d(1, sd)*x, hodge=hodge) + d(0, sd)*δ(1, sd, x, hodge=hodge)
    end

    :δ₁ => x -> inv_hodge_star(0, sd, hodge=hodge) * dual_derivative(1,sd) * ⋆(1, sd, hodge=hodge) * x
    #:i₁′ => (v,x) -> inv_hodge_star(0,sd, hodge=hodge) * wedge_product(Tuple{1,1}, sd, v, ⋆(1, sd, hodge=hodge) * x) #⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
    :i₁′ => (x,y) -> i₁′(x,y,sd,hodge)
    #:L₁′ = ... + d(0,sd)*i₁′(v,x) #i₀′(F1, d₁{X}(F1′)) + d₀{X}(i₁′(F1, F1′))
    :L₁′ => (x,y) -> L₁′(x,y,sd,hodge)
    :neg₁ => x -> -1.0 * x
    :neg₀ => x -> -1.0 * x
    :half => x -> 0.5 * x
    :third => x -> x / 3.0
    #:R₀ => x-> 1.38064852e-23 * 6.0221409e23 / (28.96 / 1000) # Boltzmann constant * ??? / (molecular mass / 1000)
    :R₀ => x-> x * boltzmann_constant * 6.0221409e23 / (mol_mass / 1000) # Boltzmann constant * ??? / (molecular mass / 1000)
    # :kᵥ => x->0.000210322*x # / density
    :kᵥ => x->1.2e-5*x # / density
    # These are the steps used to compute k.
    # We have no boundaries, so I set k to the constant k₁
    #kₜ = 0.0246295028571 # Thermal conductivity
    #k_cyl = kₜ * 4
    #density = 0.000210322
    #cₚ = 1004.703 # Specific Heat at constant pressure
    #k₁ = kₜ / (density * cₚ) # Heat diffusion constant in fluid
    #k₂ = k_cyl / (density * cₚ) # Heat diffusion constant in cylinder
    #k_col = fill(k₁, ne(s))
    #k_col[cyl] .= k₂
    #k = diagm(k_col)
    #:k => x->k₁*x
    :div₀ => (v,x) -> v / x
    :div₁ => (v,x) -> v / x
    #:avg₀₁ => x -> begin
    #  I = Vector{Int64}()
    #  J = Vector{Int64}()
    #  V = Vector{Float64}()
    #  for e in 1:ne(s)
    #      append!(J, [s[e,:∂v0],s[e,:∂v1]])
    #      append!(I, [e,e])
    #      append!(V, [0.5, 0.5])
    #  end
    #  sparse(I,J,V)*x
    #end
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(sd)
          append!(J, [sd[e,:∂v0],sd[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      sparse(I,J,V)*x
    end
    :./ => (x,y) -> x ./ y
    :.* => (x,y) -> x .* y

    # :Δ₁ => x -> begin # d ⋆ d̃ ⋆⁻¹ + ⋆ d̃ ⋆ d
    #   y = dual_derivative(0,sd)*⋆(2, sd, hodge=hodge)*d(1,sd)*x
    #   inv_hodge_star(2,sd, y; hodge=hodge) 
    #   z = d(0, sd) * inverse_hode_star(2, sd, dual_derivative(1, sd)*⋆(1,sd, hodge=hodge)*x; hodge=hodge)
    #   return y + z
    # end
    :debug => (args...)->begin println(args[1], length.(args[2:end])) end
    x=> error("Unmatched operator $my_symbol")
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length.(args))"); op(args...);end
  return (args...) ->  op(args...)
end

# End Navier Stokes

flatten_form(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

#include("coordinates.jl")
#include("spherical_meshes.jl")

radius = (6371+90) * 1e3

#primal_earth = loadmesh(ThermoIcosphere())
primal_earth = loadmesh(Icosphere(3, radius))
nploc = argmax(x -> x[3], primal_earth[:point])

orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())

physics = HeatXfer
gensim(expand_operators(physics))
sim = eval(gensim(expand_operators(physics)))

fₘ = sim(earth, generate)

begin
  vmag = 500
  # velocity(p) = vmag*ϕhat(p)
  velocity(p) = TangentBasis(CartesianPoint(p))(vmag/4, vmag/4)
  # velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, -vmag/4, 0))

# visualize the vector field
  ps = earth[:point]
  #ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)
  #GLMakie.arrows(
  #    ps, ns, fxaa=true, # turn on anti-aliasing
  #    linecolor = :gray, arrowcolor = :gray,
  #    linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
  #    align = :center, axis=(type=Axis3,)
  #)
end

#begin
v = flatten_form(velocity, earth)
#c_dist = MvNormal([radius/√(2), radius/√(2)], 20*[1, 1])
#c = 100*[pdf(c_dist, [p[1], p[2]]) for p in earth[:point]]

theta_start = 45*pi/180
phi_start = 0*pi/180
x = radius*cos(phi_start)*sin(theta_start)
y = radius*sin(phi_start)*sin(theta_start)
z = radius*cos(theta_start)
c_dist₁ = MvNormal([x, y, z], 20*[1, 1, 1])
c_dist₂ = MvNormal([x, y, -z], 20*[1, 1, 1])

c_dist = MixtureModel([c_dist₁, c_dist₂], [0.6,0.4])

t = 100*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]
t = [300.0 for p in earth[:point]] # Kelvin (300-600)
pfield = 100000*[p[3] for p in earth[:point]]

m = 8e22 # A*m^2
μ₀ = 4*π*1e-7 # kg*m/s^2/A^2
# Bθ(p) = -μ₀*m*sin(theta(p))/(4*π*radius^3)  # radius instead of r(p)*1000 
# B1Form = flatten_form(p->TangentBasis(CartesianPoint(p))(Bθ(CartesianPoint(p)), 0), earth)
# Magnetic flux density at a given point.
Br_helper(p) = -μ₀*2*m*cos(theta(p))/(4*π*radius^3)  # radius instead of r(p)*1000 
#Br_flux = hodge_star(earth, TriForm(map(triangles(earth)) do t 
#                        dual_pid = triangle_center(earth, t)
#                        p = earth[dual_pid, :dual_point]
#                        return Br(CartesianPoint(p))
#                        end))
# This is magnetic flux density.
Br_temp = inv_hodge_star(2,earth,hodge=hodge)*DualForm{0}(map(triangles(earth)) do t 
                        dual_pid = triangle_center(earth, t)
                        p = earth[dual_pid, :dual_point]
                        return Br_helper(CartesianPoint(p))
                        end)
# B₀(p) = sqrt(Bθ(p)^2+Br(p)^2)

Ii = Vector{Int64}()
Ji = Vector{Int64}()
Vi = Vector{Float64}()
for e in 1:ne(earth)
    append!(Ji, [earth[e,:∂v0],earth[e,:∂v1]])
    append!(Ii, [e,e])
    append!(Vi, [0.5, 0.5])
end
R₀ = boltzmann_constant * 6.0221409e23 / (mol_mass / 1000) # Boltzmann constant * ??? / (molecular mass / 1000)
ρ = pfield ./ (R₀*t)
momentum = v .* (sparse(Ii,Ji,Vi)*ρ)

#u₀ = construct(PhysicsState, [VectorForm(t), VectorForm(collect(v)), VectorForm(pfield)],Float64[], [:T, :V, :P])
#u₀ = construct(PhysicsState,
  #[VectorForm(Br_temp), VectorForm(t), VectorForm(collect(v)), VectorForm(pfield)],Float64[],
  #[:navierstokes_B, :continuity_T, :V, :P])
u₀ = construct(PhysicsState,
  [VectorForm(Br_temp), VectorForm(t), VectorForm(collect(v)), VectorForm(pfield), VectorForm(collect(momentum))],Float64[],
  [:navierstokes_B, :continuity_T, :V, :P, :navierstokes_M])
mesh(primal_earth, color=findnode(u₀, :P), colormap=:plasma)
tₑ = 30.0

nₑ = 1e9 # Particles per cubic meter at ~90k altitude

@info("Precompiling Solver")
my_constants = (
  navierstokes_kᵥ=1.2e-5,
  navierstokes_two=2,
  navierstokes_three=3,
  navierstokes_negone=-1,
  navierstokes_q=1.602176634e−19,
  navierstokes_nₑ=nₑ,
  navierstokes_c=299_792_458,
  continuity_diffusion_k=k₁)
prob = ODEProblem(fₘ,u₀,(0,1e-4),my_constants)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ,u₀,(0,tₑ),my_constants)
soln = solve(prob, Tsit5())
@info("Done")
#end

begin
mass(soln, t, mesh, concentration=:P) = sum(⋆(0, mesh)*findnode(soln(t), concentration))

@show extrema(mass(soln, t, earth, :P) for t in 0:tₑ/150:tₑ)
end
mesh(primal_earth, color=findnode(soln(0), :P), colormap=:jet)
mesh(primal_earth, color=findnode(soln(0) - soln(tₑ), :P), colormap=:jet)
begin
# Plot the result
times = range(0.0, tₑ, length=150)
colors = [findnode(soln(t), :P) for t in times]

# Initial frame
# fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = extrema(vcat(colors...)), colormap=:jet)
fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = (-0.0001, 0.0001), colormap=:jet)
Colorbar(fig[1,2], ob)
framerate = 5

# Animation
record(fig, "weatherNS.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = findnode(soln(t), :P)
end
end
