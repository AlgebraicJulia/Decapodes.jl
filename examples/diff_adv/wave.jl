using ACSets
using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using CairoMakie
Point2D = Point2{Float64}
Point3D = Point3{Float64}

begin
    num = 1000
    s = EmbeddedDeltaSet1D{Bool,Point2D}()
    add_vertices!(s, num, point=[Point2D(i/num * 2 * pi,0) for i in 1:num])
    add_edges!(s, 1:(nv(s)-1), 2:nv(s))
    orient!(s)
    sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point2D}(s)
    subdivide_duals!(sd, Circumcenter())
end

Wave = @decapode begin
    (U,T)::Form0
    k::Constant
    ∂ₜ(T) == -k * U
    ∂ₜ(U) == T
end

function simulate(mesh, operators, hodge = GeometricHodge())
    #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:585 =#
    #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:586 =#
    begin
        #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:174 =#
    end
    #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:587 =#
    begin
        #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:484 =#
    end
    #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:588 =#
    begin
        #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:226 =#
        var"__•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))     
        __Ṫ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
    end
    #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:589 =#
    f(du, u, p, t) = begin
            #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:589 =#
            #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:590 =#
            begin
                #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:251 =#     
                U = u.U
                T = u.T
                k = p.k
            end
            #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:591 =#
            var"•2" = Decapodes.get_tmp(var"__•2", u)
            Ṫ = Decapodes.get_tmp(__Ṫ, u)
            var"•2" .= (.-)(k)
            Ṫ .= var"•2" .* U
            #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:592 =#
            getproperty(du, :T) .= Ṫ
            getproperty(du, :U) .= T
        end
end
sim = eval(gensim(Wave))

U = map(sd[:point]) do (x,_)
    sin(x)
end

u₀ = ComponentArray(U=U,T=zeros(nv(sd)))

fₘ = simulate(sd, nothing, DiagonalHodge())

constants_and_paramters = (k = map(sd[:point]) do (x,_) x end / 3,)
tₑ = 50

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_paramters);
soln = solve(prob, Tsit5())

begin
    time = Observable(0.0)
    ys = @lift(getproperty(soln($time), :U))
    fig = lines(map(sd[:point]) do (x,_) x end, ys)
    timestamps = range(0, tₑ, step=0.2)
    record(fig, "Wave.gif", timestamps; framerate=30) do t
        time[] = t
    end
end