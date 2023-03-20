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

sim = eval(gensim(advdiffdp, [:C, :V]))
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