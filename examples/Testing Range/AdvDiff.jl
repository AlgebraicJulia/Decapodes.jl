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

function man_simulate(mesh, operators)
    begin
        #= d₀ = generate(mesh, :d₀)
        k = generate(mesh, :k)
        (⋆₁) = generate(mesh, :⋆₁)
        (-) = generate(mesh, :-)
        dual_d₁ = generate(mesh, :dual_d₁)
        (⋆₀⁻¹) = generate(mesh, :⋆₀⁻¹)=#
        L₀ = operators(mesh, :L₀)
        
        hd1 = ⋆(1,mesh,hodge=GeometricHodge())
        d0 = d(0,mesh)

        dd1 = dual_derivative(1,mesh)
        invhd0 = inv_hodge_star(0,mesh,hodge=GeometricHodge())

        ϕ₁ = Vector{Float64}(undef, ne(mesh))
        var"2" = Vector{Float64}(undef, ne(mesh))
        var"3" = Vector{Float64}(undef, ne(mesh))
        var"4" = Vector{Float64}(undef, nv(mesh))
        Ċ = Vector{Float64}(undef, nv(mesh))
        var"•1" = Vector{Float64}(undef, ne(mesh))
        ϕ₂ = Vector{Float64}(undef, ne(mesh))
        ϕ = Vector{Float64}(undef, ne(mesh))
    end
    return begin
            f(du, u, p, t) = begin
                    begin
                        C = (findnode(u, :C)).values
                        V = (findnode(u, :V)).values
                    end
                    println("--------------------")
                    # ϕ₁ = (∘(⋆₁, k, d₀))(C)
                    mul!(var"2", d0, C)
                    var"3" .= 2000 .* var"2"
                    mul!(ϕ₁, hd1, var"3")
                    var"•1" .= L₀(V, C)
                    ϕ₂ .= -var"•1"
                    ϕ .= ϕ₁ .+ ϕ₂
                    # Ċ = ((⋆₀⁻¹) ∘ dual_d₁)(ϕ)
                    mul!(var"4", dd1, ϕ)
                    mul!(Ċ, invhd0, var"4")
                    # mul!(Ċ, mat2, ϕ)
                    du .= 0.0
                    begin
                        (findnode(du, :C)).values .= Ċ
                    end
                end
        end
end

# sim = eval(gensim(advdiffdp, [:C, :V]))
# fₘ = sim(earth, generate)
fₘ = man_simulate(earth, generate)

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