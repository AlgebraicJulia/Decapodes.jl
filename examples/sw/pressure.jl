using Catlab
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

PressureFlow = quote
  # state variables
  V::Form1{X}
  P::Form0{X}
  C::Form0{X}

  # derived quantities
  ΔV::Form1{X}
  ∇P::Form1{X}
  ∇C::Form1{X}
  ΔP::Form0{X}
  ΔC::Form0{X}
  ϕₚ::Form1{X}
  ϕc::Form1{X}

  # tanvars
  V̇::Form1{X}
  Ṗ::Form0{X}
  Ċ::Form0{X}
  ∂ₜ(V) == V̇
  ∂ₜ(P) == Ṗ
  ∂ₜ(C) == Ċ
  
  ∇P == d₀(P)
  ∇C == d₀(C)
  ΔV == Δ₁(V)
  ΔP == Δ₀(P)
  ΔC == Δ₀(P)

  V̇  == α(∇P) + μ(ΔV)
  ϕₚ == γₚ(-(L₀(V, P))) 
  ϕc == γc(-(L₀(V, C))) 
  Ṗ == βₚ(Δ₀(P)) + ∘(dual_d₁,⋆₀⁻¹)(ϕₚ)
  Ċ == βc(Δ₀(C)) + ∘(dual_d₁,⋆₀⁻¹)(ϕc)

  
  # Ṗ  == ϕₚ
end

pf = SummationDecapode(parse_decapode(PressureFlow))

flatten_form(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :k => x->2000x
    :μ => x->-0.0001x
    # :μ => x->-2000x
    # :α => x->100*x
    :α => x->1.0*x  # divide alpha by length of mesh for proper scaling
    :βₚ => x->2000*x
    :γₚ => x->1*x    
    :βc => x->2000*x
    :γc => x->1*x
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
      # dual_derivative(1, sd)*⋆(1, sd)*wedge_product(Tuple{1,0}, sd, v, x)
      ⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x)
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

include("coordinates.jl")
include("spherical_meshes.jl")

radius = 6371+90

primal_earth = loadmesh(Icosphere(3,radius))
nploc = argmax(x -> x[3], primal_earth[:point])

orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3}(primal_earth)
subdivide_duals!(earth, Circumcenter())

physics = SummationDecapode(parse_decapode(PressureFlow))
gensim(expand_operators(physics), [:C, :P, :V])
sim = eval(gensim(expand_operators(physics), [:C, :P, :V]))

fₘ = sim(earth)

# begin
  vmag = 5
#   # velocity(p) = vmag*ϕhat(p)
#   # velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4, 0))
# velocity(p) = TangentBasis(CartesianPoint(p))((0, 0, 0))
velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4, 0))
#   # velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, -vmag/4, 0))

# # visualize the vector field
#   ps = earth[:point]
#   ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)
#   GLMakie.arrows(
#       ps, ns, fxaa=true, # turn on anti-aliasing
#       linecolor = :gray, arrowcolor = :gray,
#       linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
#       align = :center, axis=(type=Axis3,)
#   )
# end

begin
v = flatten_form(velocity, earth)

theta_start = 45*pi/180
phi_start = 0*pi/180
x = radius*cos(phi_start)*sin(theta_start)
y = radius*sin(phi_start)*sin(theta_start)
z = radius*cos(theta_start)
c_dist₁ = MvNormal([x, y, z], 200*[1, 1, 1])
c_dist₂ = MvNormal([x, y, -z], 200*[1, 1, 1])

c_dist = MixtureModel([c_dist₁, c_dist₂], [0.6,0.4])
# c = 100*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]
c = 10000*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]
pfield = 100000*[p[3]/radius for p in earth[:point]]
maximum(c)

u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(collect(v)), VectorForm(pfield)],Float64[], [:C, :V, :P])
fig, ax, ob = mesh(primal_earth, color=findnode(u₀, :P), colormap=:plasma)
Colorbar(fig[1,2], ob)
end
extrema(pfield)

d(0, earth, pfield)
begin
# tₑ = 30.0
tₑ = 15.0

@info("Precompiling Solver")
prob = ODEProblem(fₘ,u₀,(0,1e-4))
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")
@info("Solving")
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())
@info("Done")
end

begin
mass(soln, t, mesh, concentration=:C) = sum(⋆(0, mesh)*findnode(soln(t), concentration))

numframes = 250;

@show extrema(mass(soln, t, earth, :P) for t in 0:tₑ/numframes:tₑ)
end
mesh(primal_earth, color=findnode(soln(0), :C), colormap=:jet)
mesh(primal_earth, color=findnode(soln(0) - soln(tₑ), :C), colormap=:jet)
begin
# Plot the result
times = range(0.0, tₑ, length=numframes)
colors = [findnode(soln(t), :C) for t in times]

# Initial frame
fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = extrema(vcat(colors...)), colormap=:jet)
# fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange=(-5e-8, 5e-8), colormap=:jet)
Colorbar(fig[1,2], ob)
framerate = 5

# Animation
record(fig, "weather.gif", range(0.0, tₑ; length=numframes); framerate = 30) do t
    ob.color = findnode(soln(t), :C)
end
end

trange = range(0.0, tₑ; length=numframes)
begin
  # Plot max velocity (at any point on the grid) vs time
  maxVel = [maximum(findnode(soln(t), :V)) for t in trange]
  minVel = [minimum(findnode(soln(t), :V)) for t in trange]

end
plt = scatter(trange,maxVel)
scatter(trange,minVel,color=:red)
ylims!(minimum(minVel),maximum(maxVel))

extrema(soln(tₑ), )