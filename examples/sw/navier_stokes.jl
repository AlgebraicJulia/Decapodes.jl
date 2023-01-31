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
using GLMakie
using Logging
# using CairoMakie 

using GeometryBasics: Point3
Point3D = Point3{Float64}

# Begin Navier Stokes

# # Navier Stokes example
# DiffusionExprBody = quote
#   (T, Ṫ)::Form0{X}
#   ϕ::DualForm1{X}
#   k::Constant{X}
#   # Fick's first law
#   ϕ ==  ⋆₁(k*d₀(T))
#   # Diffusion equation
#   Ṫ == ⋆₀⁻¹(dual_d₁(ϕ))
# end

# Navier Stokes example
DiffusionExprBody = quote
  (ρ, ρ̇ )::Form0{X}
  ϕ::DualForm1{X}
  k::Constant{X}
  # Fick's first law
  ϕ ==  ⋆₁(k*d₀(ρ)) # diffusion through a dual edge
  # Diffusion equation
  ρ̇ == ⋆₀⁻¹(dual_d₁(ϕ)) # total diffusion through all dual edges about a vertex
end

Diffusion = SummationDecapode(parse_decapode(DiffusionExprBody))
to_graphviz(Diffusion)
to_graphviz(Diffusion, graph_attrs=Dict(:rankdir => "LR"))

# AdvectionExprBody = quote
#   (V)::Form1{X}  #  M = ρV
#   (ρ, p)::Form0{X}
#   Ṫ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(T))))
# end
AdvectionExprBody = quote
  (V)::Form1{X}  #  M = ρV
  (ρ)::Form0{X}
  # ρ̇ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(ρ))))
  ρ̇ == ρ*∘(⋆,d,⋆⁻¹)(V) + i₁(V, d(ρ))
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
  diffusion(ρ, ρ₁)
  advection(V, ρ, ρ₂)
  superposition(ρ, ρ̇, ρ₁, ρ₂)
end
#to_graphviz(compose_continuity, junction_labels=:variable, box_labels=:name, prog="neato")
to_graphviz(compose_continuity, junction_labels=:variable, box_labels=:name, prog="circo")

continuity_cospan = oapply(compose_continuity,
                [Open(Diffusion, [:ρ, :ρ̇ ]),
                 Open(Advection, [:V, :ρ, :ρ̇ ]),
                 Open(Superposition, [:T, :Ṫ, :Ṫ₁, :Ṫₐ])])

continuity = apex(continuity_cospan)
to_graphviz(continuity)
to_graphviz(continuity, graph_attrs=Dict(:rankdir => "LR"))

# NavierStokesExprBody = quote
#   (V, V̇, G, V)::Form1{X}
#   (T, ρ, p, ṗ)::Form0{X}
#   (two,three,kᵥ)::Constant{X}
#   V == M/avg₀₁(ρ)   # change this, also correct for mass density units
#   V̇ == neg₁(L₁′(V, V)) + 
#         kᵥ*(Δ₁(V) + d₀(δ₁(V))/three) + # include div v for compressbile fluids?
#         d₀(i₁′(V, V)/two) +
#         neg₁(d₀(p)) +    # update p with p = ρkT
#         G # + 
#         # (q*V cross product B) # 1 form in primal mesh with 0 form in dual mesh
#   ∂ₜ(V) == V̇
#   ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p)))) # *Lie(3Form) = Div(*3Form x v) --> conservation of pressure   where does this come from?
#   ∂ₜ(p) == ṗ
# end

NavierStokesExprBody = quote
  # (G, B, V)::Form1{X}
  (V)::Form1{X}
  # (ρ,  ṗ)::Form0{X}
  (two,three,kᵥ)::Constant{X}
  V̇ == neg₁(L₁′(V, V)) +     # advective term 1 of velocity
        d₀(i₁′(V, V)/two) +  # advective term 2
        kᵥ*(Δ₁(V) + d₀(δ₁(V))/three) +   # diffusive term of velocity
        neg₁(d₀(p)) 
        # G 
        # (q*V cross product B) # 1 form in primal mesh with 0 form in dual mesh
  # ∂ₜ(M) == Ṁ
  # ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p)))) # *Lie(3Form) = Div(*3Form x v) --> conservation of pressure
  # ∂ₜ(p) == ṗ
end

NavierStokes = SummationDecapode(parse_decapode(NavierStokesExprBody))
to_graphviz(NavierStokes)
to_graphviz(NavierStokes, graph_attrs=Dict(:rankdir => "LR"))

EnergyExprBody = quote
  (V)::Form1{X}
  # (p, ṗ, T, Tₑ)::Form0{X}
  (p, ṗ)::Form0{X}
  # T == p/n/kᵥ
  ṗ == 2*3 * (5/2*p*∘(⋆,d,⋆⁻¹)(V) + i₁(V, d(p))) #  - ρ*ν*3*kᵥ*(T-Tₑ) )
  
  ∂ₜ(p) == ṗ
end

Energy = SummationDecapode(parse_decapode(EnergyExprBody))
to_graphviz(Energy)
to_graphviz(Energy, graph_attrs=Dict(:rankdir => "LR"))


compose_heatXfer = @relation () begin
  continuity(V, ρ, P, T)
  navierstokes(V, ρ, P, T)
  energy(V, P)
end
#to_graphviz(compose_heatXtfer, junction_labels=:variable, box_labels=:name, prog="neato")
to_graphviz(compose_heatXfer, junction_labels=:variable, box_labels=:name, prog="circo")

heatXfer_cospan = oapply(compose_heatXfer,
                [Open(continuity, [:V, :ρ]),
                 Open(NavierStokes, [:V, :ρ, :p]),
                 Open(Energy, [:V, :p])])

HeatXfer = apex(heatXfer_cospan)
to_graphviz(HeatXfer)
to_graphviz(HeatXfer, graph_attrs=Dict(:rankdir => "LR"))

hodge = GeometricHodge()
d₀(x,sd) = d(0,sd)*x
d₁(x,sd) = d(1,sd)*x
⋆₀(x,sd,hodge) = ⋆(0,sd,hodge=hodge)*x
⋆₁(x,sd,hodge) = ⋆(1,sd,hodge=hodge)*x
⋆₂(x,sd,hodge) = ⋆(2,sd,hodge=hodge)*x
⋆₀⁻¹(x,sd,hodge) = inv_hodge_star(0,sd,hodge=hodge)*x
⋆₁⁻¹(x,sd,hodge) = inv_hodge_star(1,sd,hodge=hodge)*x
∧₀₁′(x,y,sd) = wedge_product(Tuple{0,1}, sd, x, y)
∧₁₁′(x,y,sd) = wedge_product(Tuple{1,1}, sd, x, y)
#∧₁₀′(x,y,sd) = cp_2_i ...
i₀′(x,y,sd,hodge) = -1.0 * (∧₁₀′(x, ⋆₂(y,sd,hodge)),sd)
#i₁′(v,x) = inv_hodge_star(0,sd, hodge=hodge) * wedge_product(Tuple{1,1}, sd, v, ⋆(1, sd, hodge=hodge) * x) #⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
i₁′(x,y,sd,hodge) = ⋆₀⁻¹(∧₁₁′(x, ⋆₁(y,sd,hodge), sd),hodge) #⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
L1′(x,y,sd,hodge) = i₀′(x,d₁(y,sd),sd,hodge) + d₀(i₁′(x,y,sd,hodge),sd)
boltzmann_constant = 1.38064852e-23
mol_mass = 28.96
density = 0.000210322
kₜ = 0.0246295028571 # Thermal conductivity
cₚ = 1004.703 # Specific Heat at constant pressure
k₁ = kₜ / (density * cₚ) # Heat diffusion constant in fluid

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :k => x->2000x
    :μ => x->-0.0001x
    # :μ => x->-2000x
    :β => x->2000*x
    :γ => x->1*x
    :⋆₀ => x->⋆(0,sd,hodge=hodge)*x
    :⋆₁ => x->⋆(1, sd, hodge=hodge)*x
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=hodge)
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=hodge)*x
    :d₀ => x->d(0,sd)*x
    :d₁ => x->d(1,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
    :plus => (+)
    :(-) => x-> -x
    # :L₀ => (v,x)->dual_derivative(1, sd)*(i0(v, x))
    :L₀ => (v,x)->begin
      dual_derivative(1, sd)*⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x) # wedge product of 1 and 2 form???  switched with 183; add in lv*d term
      # ⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x)
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
    :R₀ => x-> boltzmann_constant * 6.0221409e23 / (mol_mass / 1000) # Boltzmann constant * ??? / (molecular mass / 1000)
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
    :k => x->k₁*x
    :div₀ => (v,x) -> v / x
    :div₁ => (v,x) -> v / x
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(s)
          append!(J, [s[e,:∂v0],s[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      sparse(I,J,V)*x
    end

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

radius = (6371+90) * 1e3  # to shell altitude: meters

#primal_earth = loadmesh(ThermoIcosphere())
primal_earth = loadmesh(Icosphere(3, radius))
nploc = argmax(x -> x[3], primal_earth[:point])

orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())

#physics = SummationDecapode(parse_decapode(PressureFlow))
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
  ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)
  GLMakie.arrows(
      ps, ns, fxaa=true, # turn on anti-aliasing
      linecolor = :gray, arrowcolor = :gray,
      linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
      align = :center, axis=(type=Axis3,)
  )
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
pfield = 100000*[p[3] for p in earth[:point]]

m = 8e22 # dipole moment units: A*m^2
μ₀ = 4*π*1e-7 # permeability units: kg*m/s^2/A^2
# Bθ(p) = -μ₀*m*sin(theta(p))/(4*π*radius^3)  # radius instead of r(p)*1000 
# B1Form = flatten_form(p->TangentBasis(CartesianPoint(p))(Bθ(CartesianPoint(p)), 0), earth)
Br(p) = -μ₀*2*m*cos(theta(p))/(4*π*radius^3)  # radius instead of r(p)*1000 
Br_flux = hodge_star(earth, TriForm(map(triangles(earth)) do t 
                        dual_pid = triangle_center(earth, t)
                        p = earth[dual_pid, :dual_point]
                        return Br(CartesianPoint(p))
                        end))
# B₀(p) = sqrt(Bθ(p)^2+Br(p)^2)


u₀ = construct(PhysicsState, [VectorForm(t), VectorForm(collect(v)), VectorForm(pfield)],Float64[], [:T, :V, :P])
mesh(primal_earth, color=findnode(u₀, :P), colormap=:plasma)
tₑ = 30.0

@info("Precompiling Solver")
my_constants = (kᵥ=1.2e-5,two=2,three=3,continuity_diffusion_k=k₁)
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

#AdvDiff = quote
#    C::Form0{X}
#    Ċ::Form0{X}
#    V::Form1{X}
#    ϕ::Form1{X}
#    ϕ₁::Form1{X}
#    ϕ₂::Form1{X}
#    starC::DualForm2{X}
#    lvc::Form1{X}
#    # Fick's first law
#    ϕ₁ ==  ∘(d₀,k,⋆₁)(C)
#    # Advective Flux
#    ϕ₂ == -(L₀(V, C))
#    # Superposition Principle
#    ϕ == plus(ϕ₁ , ϕ₂)
#    # Conservation of Mass
#    Ċ == ∘(dual_d₁,⋆₀⁻¹)(ϕ)
#    ∂ₜ(C) == Ċ
#end
#
#NavierStokes = quote
#  V::Form1{X}
#  V̇::Form1{X}
#  G::Form1{X}
#  T::Form0{X}
#  ρ::Form0{X}
#  ṗ::Form0{X}
#  p::Form0{X}
#  
#
#  V̇ == neg₁(L₁′(V, V)) + 
#        div₁(kᵥ(Δ₁(V) + third(d₀(δ₁(V)))), avg₀₁(ρ)) +
#        d₀(half(i₁′(V, V))) +
#        neg₁(div₁(d₀(p),avg₀₁(ρ))) +
#        G
#  ∂ₜ(V) == V̇
#  ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
#  ∂ₜ(p) == ṗ
#end
#
#parse_decapode(NavierStokes)
#SummationDecapode(parse_decapode(NavierStokes))

# Energy = @decapode Flow2DQuantities begin
#   (V)::Form1{X}
#   (ρ, p, T, Ṫ, Ṫₐ, Ṫ₁, bc₀)::Form0{X}

#   ρ == div₀(p, R₀(T))
#   Ṫₐ == neg₀(⋆₀⁻¹{X}(L₀(V, ⋆₀{X}(T))))
#   ∂ₜₐ(Ṫₐ) == bc₀
#   Ṫ == Ṫₐ + Ṫ₁
#   ∂ₜ{Form0{X}}(T) == Ṫ
# end

# BoundaryConditions = @decapode Flow2DQuantities begin 
#   (V, V̇, bc₁)::Form1{X}
#   (Ṫ, ṗ, bc₀)::Form0{X} 
#   # no-slip edges
#   ∂ᵥ(V̇) == bc₁
#   # No change on left/right boundaries
#   ∂ₜ(Ṫ) == bc₀
#   ∂ₚ(ṗ) == bc₀
# end
