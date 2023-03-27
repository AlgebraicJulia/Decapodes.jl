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
using TerminalLoggers
using Logging
using CUDA
using StaticArrays
global_logger(TerminalLogger())

Point3D = Point3{Float64}

flatten(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
        _ => default_dec_generate(sd, my_symbol, hodge)
    end
  
    return (args...) ->  op(args...)
end

begin
    primal_earth = loadmesh(Icosphere(5))
    primal_earth = EmbeddedDeltaSet2D("Icosphere8.obj")
    nploc = argmax(x -> x[3], primal_earth[:point])
    primal_earth[:edge_orientation] .= false
    orient!(primal_earth)
    earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth)
    subdivide_duals!(earth, Circumcenter());
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
  
function man_simulate(mesh, operators)
    begin
        # Δ₀ = generate(mesh, :Δ₀)
        mesh = earth
        lpdr0   = Δ(0, mesh)
        var"•2" = Vector{Float64}(undef, nv(mesh))
        U2V     = Vector{Float64}(undef, nv(mesh))
        var"•1" = Vector{Float64}(undef, nv(mesh))
        aTU     = Vector{Float64}(undef, nv(mesh))
        var"•4" = Vector{Float64}(undef, nv(mesh))
        var"•6" = Vector{Float64}(undef, nv(mesh))
        var"•5" = Vector{Float64}(undef, nv(mesh))
        sum_1   = Vector{Float64}(undef, nv(mesh))
        var"•3" = Vector{Float64}(undef, nv(mesh))
        #U̇ = Vector{Float64}(undef, nv(mesh))
        #V̇ = Vector{Float64}(undef, nv(mesh))
    end
    return begin
        f(du, u, p, t) = begin
            begin
                U = (findnode(u, :U)).values
                V = (findnode(u, :V)).values
                One = p.One
                fourfour = p.fourfour
                threefour = p.threefour
                α = p.α
                F = p.F(t)
            end
            # println("--------------------")
            # var"•2" = Δ₀(U)
            mul!(var"•2", lpdr0, U)
            var"•1" .= U .* U
            U2V .= var"•1" .* V
            aTU .= α .* var"•2"
            var"•4" .= fourfour .* U
            var"•6" .= threefour .* U
            var"•5" .= var"•6" .- U2V
            sum_1 .= One .+ U2V
            (findnode(du, :V)).values .= var"•5" .+ aTU
            var"•3" .= sum_1 .- var"•4"
            (findnode(du, :U)).values .= var"•3" .+ aTU .+ F
        end
    end
end

# sim = eval(gensim(Brusselator))
# fₘ = sim(earth, generate)
fₘ = man_simulate(earth, generate)

begin
    U = map(earth[:point]) do (_,y,_)
      abs(y)
    end
    
    V = map(earth[:point]) do (x,_,_)
      abs(x)
    end

    One = ones(nv(earth))
        
    F₁ = map(earth[:point]) do (_,_,z)
      z ≥ 0.8 ? 5.0 : 0.0
    end

    F₂ = zeros(nv(earth))

    constants_and_parameters = (
        fourfour = 4.4,
        threefour = 3.4,
        α = 0.001,
        F = t -> t ≥ 1.1 ? F₁ : F₂,
        One = One)

    u₀ = construct(PhysicsState, [VectorForm(U), VectorForm(V)],Float64[], [:U, :V])
    tₑ = 11.5
    prob = ODEProblem(fₘ,u₀,(0, tₑ), constants_and_parameters)
    soln = solve(prob, Tsit5(), save_everystep = false)
end

fig, ax, ob = GLMakie.mesh(primal_earth, color = findnode(soln(tₑ), :U))
for t in range(0.0, tₑ; length=300)
    sleep(0.1)
    ob.color = findnode(soln(t), :U)
end


#UnicodePlots.scatterplot(soln.t[begin+1:end] - soln.t[begin:end-1])
     ┌─────────────────────────────────────────────┐ 
0.07 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡏⠂⠀⠀⠀⠀⠀⠀⠀⠀⢀⡄⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⢡⣄⣠⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⡤⠤⠖⠁⣦⣀⣀⣀⣀⣀⣀⣀⣀⡼⣁⣀⣀⡀⠀⠀⠀⠀⠀│ 
     │⠚⠀⢉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠚⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀│ 
   0 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
     └─────────────────────────────────────────────┘ 
     ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀600⠀ 


# GPU on Icosphere 5 Brusselator for 15.0 sim time.
#     BenchmarkTools.Trial: 8 samples with 1 evaluation.
# Range (min … max):  622.786 ms … 668.150 ms  ┊ GC (min … max): 0.00% … 0.00%
# Time  (median):     633.515 ms               ┊ GC (median):    0.00%
# Time  (mean ± σ):   639.551 ms ±  15.785 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#  █         ██ █ █   █                             █          █
#  █▁▁▁▁▁▁▁▁▁██▁█▁█▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁█ ▁
#  623 ms           Histogram: frequency by time          668 ms <
#
# Memory estimate: 22.76 MiB, allocs estimate: 391054.

# GPU on Icosphere 5 Brusselator for 15.0 sim time, w/ @inbounds
#BenchmarkTools.Trial: 8 samples with 1 evaluation.
# Range (min … max):  603.889 ms … 655.395 ms  ┊ GC (min … max): 0.00% 
#… 0.00%
# Time  (median):     633.734 ms               ┊ GC (median):    0.00% 
# Time  (mean ± σ):   633.116 ms ±  15.405 ms  ┊ GC (mean ± σ):  0.00% 
#± 0.00%
#
#  █                     █       █  █  █       █  █            █       
#  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁█▁▁█▁▁▁▁▁▁▁█▁▁█▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
#  604 ms           Histogram: frequency by time          655 ms <     
#
# Memory estimate: 22.76 MiB, allocs estimate: 391115.

# CPU on Icosphere 5 Brusselator for 15.0 sim time.
#BenchmarkTools.Trial: 45 samples with 1 evaluation.
# Range (min … max):   65.203 ms … 195.998 ms  ┊ GC (min … max):  0.00% … 47.45%
# Time  (median):     108.710 ms               ┊ GC (median):     0.00% Time  (mean ± σ):   111.422 ms ±  30.407 ms  ┊ GC (mean ± σ):  12.80% ± 18.05%
#
#  ▄                  ▄▄█ 
#  █▆▁▄▁▁▁▁▁▁▁▁▁▄▁▁▆█▆███▄▁▁▁▄▁▄▁▁▁▁▁▄▄▄▁▁▁▁▄▆▁▁▁▁▁▄▁▁▁▁▄▁▁▁▁▁▁▄ ▁     
#  65.2 ms          Histogram: frequency by time          196 ms <     
#
# Memory estimate: 70.99 MiB, allocs estimate: 12102.