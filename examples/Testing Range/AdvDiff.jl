using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
using GeometryBasics: Point3

Point3D = Point3{Float64}

flatten(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
      :k => x->2000x
      _ => default_dec_generate(sd, my_symbol, hodge)
    end
  
    return (args...) ->  op(args...)
  end
  
begin
    RADIUS = 6371+90
    primal_earth = loadmesh(Icosphere(4, RADIUS))
    nploc = argmax(x -> x[3], primal_earth[:point])
    primal_earth[:edge_orientation] .= false
    orient!(primal_earth)
    earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
    subdivide_duals!(earth, Circumcenter())
end

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
    advdiffdp = SummationDecapode(advdiff)
end

#function man_simulate(mesh, operators)
#    begin
#        d₀ = generate(mesh, :d₀)
#        (⋆₁) = generate(mesh, :⋆₁)
#        dual_d₁ = generate(mesh, :dual_d₁)
#        (⋆₀⁻¹) = generate(mesh, :⋆₀⁻¹)
#        L₀ = operators(mesh, :L₀)
#        k = generate(mesh, :k)
#        (-) = generate(mesh, :-)
#        plus = operators(mesh, :plus)
#        
#        # TODO: Could place matrix allocation code in here instead of in generate
#        # Would allow for composite functions to use previously found matrices
#
#        #= tmpstar1 = ⋆(1, earth, hodge=GeometricHodge())
#        tmpinvstar1 = inv_hodge_star(0,earth; hodge=GeometricHodge())
#        tmpd0 = d(0,earth)
#        tmpdd1 = dual_derivative(1,earth)
#
#        ⋆₁ = x -> tmpstar1 * x
#        ⋆₀⁻¹ = x -> tmpinvstar1 * x
#        d₀ = x -> tmpd0 * x
#        dual_d₁ = x -> tmpdd1 * x
#        L₀ = (v, x) -> tmpstar1*wedge_product(Tuple{1,0}, earth, v, x) =#
#
#        tmpstar1 = ⋆(1, mesh, hodge=GeometricHodge())
#        val_pack = p_wedge_product(1, mesh)
#        L₀ = (v, x) -> tmpstar1 * c_wedge_product(1, x, v, val_pack)
#
#        end
#    return begin
#        f(du, u, p, t) = begin
#            begin
#                C = (findnode(u, :C)).values
#                V = (findnode(u, :V)).values
#            end
#            var"•_1_1" = d₀(C)
#            var"•_1_2" = k(var"•_1_1")
#            ϕ₁ = ⋆₁(var"•_1_2")
#            var"•1" = L₀(V, C)
#            ϕ₂ = -var"•1"
#            ϕ = plus(ϕ₁, ϕ₂)
#            var"•_3_1" = dual_d₁(ϕ)
#            Ċ = ⋆₀⁻¹(var"•_3_1")
#            du .= 0.0
#            begin
#                (findnode(du, :C)).values .= Ċ
#            end
#        end
#    end
#end

# TODO: Code actually just uses this simulate function, not the above
#fₘ = eval(man_simulate(earth, generate))
 sim = eval(gensim(AdvDiff, [:C, :V]))
 fₘ = sim(earth, generate)

 begin
    c_dist = MvNormal(nploc[[1,2]], 100[1, 1])
    c = [pdf(c_dist, [p[1], p[2]]./√RADIUS) for p in earth[:point]]

    vmag = 500
    velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, vmag/4))
    v = flatten(velocity, earth)

    u₀ = construct(PhysicsState, [VectorForm(c), VectorForm(collect(v))],Float64[], [:C, :V])
    tₑ = 10
    prob = ODEProblem(fₘ,u₀,(0, tₑ))
    soln = solve(prob, Tsit5())
end

fig, ax, ob = GLMakie.mesh(primal_earth, color = findnode(soln(0), :U))
for t in range(0.0, tₑ; length=150)
  sleep(0.05)
  ob.color = findnode(soln(t), :U)
end