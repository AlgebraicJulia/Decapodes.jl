using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
using ComponentArrays
using GeometryBasics: Point3

Point3D = Point3{Float64}

""" Wedge product of a 0-form and a ``k``-form.
"""
function my_wedge_product_zero(::Type{Val{k}}, s::HasDeltaSet,
                            f, α, x::Int) where k
  subs = subsimplices(k, s, x)
  vs = primal_vertex(k, s, subs)
  coeffs = map(x′ -> dual_volume(k,s,x′), subs) / volume(k,s,x)
  dot(coeffs, f[vs]) * α[x] / factorial(k)
end

my_wp(::Type{Tuple{0,k}}, s::HasDeltaSet, f, β, x::Int) where k =
  my_wedge_product_zero(Val{k}, s, f, β, x)


""" Alias for the wedge product operator [`∧`](@ref).
"""
const my_wedge_product = my_wp

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
    :∧₀₁ => (x,y)-> map(simplices(0+1, sd)) do z
                        ##∧(Tuple{0,1}, sd, x, y, z)

                        subs = subsimplices(1, sd, z)
                        vs = primal_vertex(1, sd, subs)
                        coeffs = map(x′ -> dual_volume(1,sd,x′), subs) / volume(1,sd,z)
                        dot(coeffs, x[vs]) * y[z] / factorial(1)
                    end
    :plus => (+)
  end
  ## return (args...) -> begin println("applying $my_symbol"); println("arg length $(length(args[1]))"); op(args...);end
  return (args...) ->  op(args...)
end


AdvectionDiffusionExprBody =  quote
    C::Form0{X}
    Ċ::Form0{X}
    V::Form1{X}
    ϕ::Form1{X}
    ϕ₁::Form1{X}
    ϕ₂::Form1{X}

    ## Fick's first law
    ϕ₁ ==  (d₀∘k)(C)
    ϕ₂ == ∧₀₁(C,V)
    ϕ == ϕ₁ + ϕ₂
    ## Diffusion equation
    Ċ == ∘(⋆₁, dual_d₁, ⋆₀⁻¹)(ϕ)
    ∂ₜ(C) == Ċ
end


diffExpr = parse_decapode(AdvectionDiffusionExprBody)
ddp = SummationDecapode(diffExpr)
gensim(expand_operators(ddp), [:C, :V])
f = eval(gensim(expand_operators(ddp), [:C, :V]))

radius = 6371+90
primal_earth = loadmesh(Icosphere(4, radius))
nploc = argmax(x -> x[3], primal_earth[:point])
orient!(primal_earth)
earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
subdivide_duals!(earth, Circumcenter())

fₘ = f(earth, generate)
c_dist = MvNormal(nploc[[1,2]], 100[1, 1])
c = [pdf(c_dist, [p[1], p[2]]./√radius) for p in earth[:point]]
v = ones(Float64, ne(earth))

wedge_product(Tuple{0,1}, earth, c, v)

u₀ = ComponentArrays(C=c, V=v)
tₑ = 1000
prob = ODEProblem(fₘ,u₀,(0,tₑ))
soln = solve(prob, Tsit5())

using CairoMakie

mesh(primal_earth, color=soln(0).C, colormap=:plasma)
mesh(primal_earth, color=soln(tₑ).C, colormap=:plasma)
mesh(primal_earth, color=soln(tₑ)-soln(0).C, colormap=:plasma)

