using Catlab
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form
using ComponentArrays
using DiagrammaticEquations
using DiagrammaticEquations.Deca
using Decapodes
using StaticArrays
using ComponentArrays
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis
#
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}
Stokes_Chorin = @decapode begin
    # P := Press𝐮re.
    # Pₙ := Press𝐮re at the next time-step.
    (P,Pₙ)::Form0
    # 𝐮 := Velocity.
    # (Note: We assume ρ is 1.)
    𝐮::Form1
    # μ := Viscosity.
    # Δt := Time-step.
    (φ,μ,Δt)::Constant
    𝐮ᵢ == 𝐮 + (μ * Δ(𝐮))*Δt
    Pₙ == Δ⁻¹(∘(⋆,d,⋆)(𝐮 / Δt))
    ∂ₜ(𝐮) == (𝐮ᵢ - Δt*d(Pₙ) - 𝐮)/Δt
    ∂ₜ(P) == (Pₙ - P)/Δt
end; infer_types!(Stokes_Chorin)
s = triangulated_grid(1,1,1/100,1/100,Point3D);
s[:point] = map(s[:point]) do (x,y,z)
    Point3D(-y,x,z)
end
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());
# functions

function generate(sd, my_symbol; hodge=DiagonalHodge())
    Δ0 = Δ(0,sd);
    fΔ0 = LinearAlgebra.factorize(Δ0);
    op = @match my_symbol begin
        :Δ⁻¹ => x -> begin
            y = fΔ0 \ x
            y .- minimum(y)
        end
        :.* => (x,y) -> x .* y
        _ => error("$my_symbol is not a supported operator.")
    end
    return (args...) -> op(args...)
end

symsim = gensim(Stokes_Chorin);
open("Stokes_Chorin.jl", "w") do f
    write(f, string(symsim))
end
sim = include("Stokes_Chorin.jl")
f = sim(sd, generate);

Δt=1e-5
ω = eval_constant_primal_form(sd, SVector{3}([1e-6,1e-6,0.]));
u₀ = ComponentArray(𝐮=ω,P=2*ones(nv(sd)))
du = similar(u₀);
constants_and_parameters=ComponentArray(Δt=Δt,μ=1,φ=[1.0 for _ in 1:length(ω)],lbP=lbP, rbP=rbP, tbP=tbP, bbP=bbP)

#t=1e-5;
t=1e-6;
prob = ODEProblem(f,u₀,(0,t),constants_and_parameters);
soln = solve(prob, Tsit5(), adaptive=false, dt=Δt);
soln.t

function juxtapose(pode, msh, mshd, soln, vars; time::Float64=0)
    sm = ♯_mat(mshd, AltPPSharp())
    fig = Figure()
    foreach(enumerate(vars)) do (k,var)
        ax = CairoMakie.Axis(fig[1,2k-1], title=string(var))
        f = @match only(pode[incident(pode, var, :name), :type]) begin
            :Form1 => u -> norm.(sm * u)
            _ => identity
        end
        _m=mesh!(ax, msh, color=f(getproperty(soln(time), var)), colormap=:jet)
        Colorbar(fig[1,2k], _m)
    end
    fig
end

juxtapose(Stokes_Chorin, s, sd, soln, [:P, :𝐮], time=t)
