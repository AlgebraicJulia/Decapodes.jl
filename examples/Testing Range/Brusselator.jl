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
    primal_earth = loadmesh(Icosphere(2, RADIUS))
    nploc = argmax(x -> x[3], primal_earth[:point])
    primal_earth[:edge_orientation] .= false
    orient!(primal_earth)
    earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
    subdivide_duals!(earth, Circumcenter())
end

begin
    Brusselator = SummationDecapode(parse_decapode(quote
        (U, V)::Form0{X} 
        (U2V, One, aTU)::Form0{X}
        (U̇, V̇)::Form0{X}

        (fourfour, threefour, α)::Constant{X}
        F::Parameter{X}

        U2V == (U .* U) .* V
        aTU == α * Δ(U)
        
        U̇ == One + U2V - (fourfour * U) + aTU + F
        V̇ == (threefour * U) - U2V + aTU

        ∂ₜ(U) == U̇
        ∂ₜ(V) == V̇
    end))

    bespoke_op1_inf_rules = [
        (src_type = :Form0, tgt_type = :infer, replacement_type = :Form0, op = :Δ)]
  
    bespoke_op2_inf_rules = [
        (proj1_type = :Form0, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :.*),
        (proj1_type = :Form0, proj2_type = :Parameter, res_type = :infer, replacement_type = :Form0, op = :*),
        (proj1_type = :Parameter, proj2_type = :Form0, res_type = :infer, replacement_type = :Form0, op = :*)]
  
    infer_types!(Brusselator,
        vcat(bespoke_op1_inf_rules, op1_inf_rules_2D),
        vcat(bespoke_op2_inf_rules, op2_inf_rules_2D))

    resolve_overloads!(Brusselator)
end
  
sim = eval(gensim(Brusselator))
fₘ = sim(earth, generate)
 
begin
    sd = earth
    U = map(sd[:point]) do (_,y,_)
      abs(y)
    end
    
    V = map(sd[:point]) do (x,_,_)
      abs(x)
    end
    
    # TODO: Try making this sparse.
    F₁ = map(sd[:point]) do (_,_,z)
      z ≥ 0.8 ? 5.0 : 0.0
    end

    F₂ = zeros(nv(sd))

    One = ones(nv(sd))

    constants_and_parameters = (
        fourfour = 4.4,
        threefour = 3.4,
        α = 0.01,
        F = t -> t ≥ 1.1 ? F₁ : F₂)

    u₀ = construct(PhysicsState, [VectorForm(U), VectorForm(V), VectorForm(One)],Float64[], [:U, :V, :One])
    tₑ = 0.01
    prob = ODEProblem(fₘ,u₀,(0, tₑ), constants_and_parameters)
    soln = solve(prob, Tsit5())
end
