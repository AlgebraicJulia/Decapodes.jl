using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra

function generate(sd, my_symbol; hodge=DiagonalHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :k => x->2000x
    :⋆₀ => x->⋆(0,sd,hodge=hodge)*x
    :⋆₁ => x->⋆(1, sd, hodge=hodge)*x
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=hodge)
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=hodge)*x
    :d₀ => x->d(0,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
    :plus => (+)
    :(-) => x-> -x
    # :L₀ => (v,x)->dual_derivative(1, sd)*(i0(v, x))
    :L₀ => (v,x)->begin
      # dual_derivative(1, sd)*⋆(1, sd)*wedge_product(Tuple{1,0}, sd, v, x)
      ⋆(1, sd)*wedge_product(Tuple{1,0}, sd, v, x)
    end
    :i₀ => i0 
    :debug => (args...)->begin println(args[1], length.(args[2:end])) end
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length.(args))"); op(args...);end
  return (args...) ->  op(args...)
end


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

include("spherical_meshes.jl")
radius = 6371+90
primal_earth, npi, spi = makeSphere(0, 180, 5, 0, 360, 5, radius);
nploc = primal_earth[npi, :point]
orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())


fₘ = f(earth)
c_dist = MvNormal(nploc[[1,2]], 100[1, 1])
c = [pdf(c_dist, [p[1], p[2]]./√radius) for p in earth[:point]]

u₀ = construct(PhysicsState, [VectorForm(c)],Float64[], [:C])
tₑ = 10
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

# using CairoMakie 
using GLMakie

mesh(primal_earth, color=findnode(soln(0), :C), colormap=:plasma)
mesh(primal_earth, color=findnode(soln(tₑ), :C), colormap=:plasma)
mesh(primal_earth, color=findnode(soln(tₑ)-soln(0), :C), colormap=:plasma)

begin
AdvDiff = quote
    C::Form0{X}
    Ċ::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}
    starC::DualForm2{X}
    lvc::Form1{X}
    # Fick's first law
    ϕ₁ ==  ∘(d₀,k,⋆₁)(C)
    # Advective Flux
    ϕ₂ == -(L₀(V, C))
    # Superposition Principle
    ϕ == plus(ϕ₁ , ϕ₂)
    # Conservation of Mass
    Ċ == ∘(dual_d₁,⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end

advdiff = parse_decapode(AdvDiff)
advdiffdp = NamedDecapode(advdiff)
gensim(expand_operators(advdiffdp), [:C, :V])
end 
begin
sim = eval(gensim(expand_operators(advdiffdp), [:C, :V]))

fₘ = sim(earth)

# velocity(p) = [-p[2]/p[1], 1.0, 0]/log(abs(p[3])+1)
phi(p) = atan(p[2]/p[1])
theta(p) = atan(sqrt(p[2]^2 + p[1]^2)/p[3])
θhat(p) = [cos(phi(p))*cos(theta(p)), sin(phi(p))*cos(theta(p)), -sin(theta(p))]
ϕhat(p) = [-sin(phi(p)), cos(phi(p)), 0]
rhat(p) = [cos(phi(p))*sin(theta(p)), sin(phi(p))*sin(theta(p)), cos(theta(p))]
v = -1
velocity(p) = v*ϕhat(p)#/log(abs(p[3])+1)
# velocity(p) = [-p[2]/p[1], 0.0, sign(p[1]*abs(p[3]))]#/log(abs(p[3])+1)
# velocity(p) = begin
#   θ = atan(-p[2]/p[1])
#   return [sin(θ), cos(θ), 0]*abs(norm(p)-p[3])/norm(p)
# end
# velocity(p) = [0, 0, sign(p[1]*abs(p[3]))]#/log(abs(p[3])+1)
v = ♭(earth, DualVectorField(velocity.(earth[triangle_center(earth),:dual_point])))
theta_start = 45*pi/180
phi_start = 0*pi/180
x = radius*cos(phi_start)*sin(theta_start)
y = radius*sin(phi_start)*sin(theta_start)
z = radius*cos(theta_start)
c_dist = MvNormal([x, y, z], 20*[1, 1, 1])
c = 100*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]

u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(5000v)],Float64[], [:C, :V])
mesh(primal_earth, color=findnode(u₀, :C), colormap=:jet)
tₑ = 2.2
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())
end

Figure()
scatter(0:tₑ/150:tₑ, [norm(findnode(soln(t), :C)) for t in 0:tₑ/150:tₑ ])

mesh(primal_earth, color=findnode(soln(tₑ), :C), colormap=:jet)
begin
# Plot the result
times = range(0.0, tₑ, length=150)
colors = [findnode(soln(t), :C) for t in times]

# Initial frame
# fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = extrema(vcat(colors...)), colormap=:jet)
fig, ax, ob = mesh(primal_earth, color=colors[1], colorrange = (-0.0001, 0.0001), colormap=:jet)
Colorbar(fig[1,2], ob)
framerate = 5

# Animation
record(fig, "diff_adv.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    ob.color = findnode(soln(t), :C)
end
end

cords(mesh) = begin
  dim(mesh, i) = [p[i] for p in unique(mesh[:point]) if p[3] > 0]
  return dim.([mesh], [1,2,3])
end
# surface(cords(earth)..., axis=(type=Axis3,))
# scene = Scene();
# arr = Makie.arrows!(
#     cords(earth)..., ones(nv(earth)), ones(nv(earth)), ones(nv(earth));
#     arrowsize = 10.1, linecolor = (:gray, 0.7), linewidth = 0.02, lengthscale = 0.1
# )

# velocity(p) = [-p[2]/p[1], 1.0, sign(p[1]*abs(p[3]))]#/log(abs(p[3])+1)
begin
  ϕ(p) = atan(p[2]/p[1])
  θ(p) = atan(sqrt(p[2]^2 + p[1]^2)/p[3])
  r(p) = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
  θhat(p) = [cos(ϕ(p))*cos(θ(p)), sin(ϕ(p))*cos(θ(p)), -sin(θ(p))]
  ϕhat(p) = [-sin(ϕ(p)), cos(ϕ(p)), 0]
  rhat(p) = [cos(ϕ(p))*sin(θ(p)), sin(ϕ(p))*sin(θ(p)), cos(θ(p))]
  v = -1
  velocity(p) = [v*θhat(p)]#/log(abs(p[3])+1)
  
ps = earth[:point]
# ns = 100*Vec3f.(map(velocity, ps))
arrows(
    ps, ns, fxaa=true, # turn on anti-aliasing
    linecolor = :gray, arrowcolor = :gray,
    linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
    align = :center, axis=(type=Axis3,)
)
end