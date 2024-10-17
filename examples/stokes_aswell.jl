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

StokesDynamics = @decapode begin
    P::Form0 ## Pressure.
    v::Form1 ## Velocity.
    (φ,μ)::Constant
    # Associate tangent variables with a state variable.
    ∂ₜ(v) == v̇
    ∂ₜ(P) == Ṗ
    # the equations
    v̇ == μ * Δ(v)-d₀(P) + φ
    Ṗ == ⋆₀⁻¹(dual_d₁(⋆₁(v)))
end; infer_types!(StokesDynamics)
# I wonder if there is a principled way of creating the boundary pode
StokesBoundaries = @decapode begin
    lbP::Constant
    rbP::Constant
    tbP::Constant
    bbP::Constant
end
StokesMorphism = @relation () begin
    LeftBoundary(C, lbC)
    RightBoundary(C, rbC)
    TopBoundary(C, tbC)
    BottomBoundary(C, bbC)
end
StokesCollage = DiagrammaticEquations.collate(
    StokesDynamics,
    StokesBoundaries,
    StokesMorphism,
    Dict(:C => :Ṗ,
         :lbC => :lbP,
         :rbC => :rbP,
         :tbC => :tbP,
         :bbC => :bbP))

# mesh
s = triangulated_grid(1,1,1/100,1/100,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());

# functions
diam(s) = sqrt(sum((s[:point][s[:∂v0][1]]- s[:point][s[:∂v1][1]]) .^2))
left_wall_idxs(s) = findall(x -> abs(x[1]) < diam(s), s[:point])
right_wall_idxs(s) = findall(x -> abs(x[1]-1) < diam(s), s[:point])
bottom_wall_idxs(s) = findall(x -> abs(x[2]) < diam(s), s[:point])
top_wall_idxs(s) = findall(x -> abs(x[2]-1) < diam(s), s[:point])
function generate(sd, my_symbol; hodge=GeometricHodge())
    op = @match my_symbol begin
        :LeftBoundary => (ℓ, ℓ₀) -> begin
            ℓ[left_wall_idxs(sd)] .= ℓ₀
            ℓ
        end
        :RightBoundary => (ℓ, ℓ₀) -> begin
            ℓ[right_wall_idxs(sd)] .= ℓ₀
            ℓ
        end
        :TopBoundary => (ℓ, ℓ₀) -> begin
            ℓ[top_wall_idxs(sd)] .= ℓ₀
            ℓ
        end
        :BottomBoundary => (ℓ, ℓ₀) -> begin
            ℓ[bottom_wall_idxs(sd)] .= ℓ₀
            ℓ
        end
        # :Boundary => (ℓ,ℓ₀) -> begin
            # ℓ[left_wall_idxs(sd)] .= ℓ₀
            # ℓ[right_wall_idxs(sd)] .= ℓ₀
        :.* => (x,y) -> x .* y
        _ => error("A!!!!!!!!!")
    end
    return (args...) -> op(args...)
end

init() = ComponentArray(idxs=Vector{Int}(), vals=Vector{Float64}())
ComponentArray(l=init(),r=init())

lbP = fill(1.0, length(left_wall_idxs(sd)));
rbP = fill(1.0, length(right_wall_idxs(sd)));
tbP = fill(1.0, length(top_wall_idxs(sd)));
bbP = fill(1.0, length(bottom_wall_idxs(sd)));

# initialize data
symsim = gensim(StokesCollage);
open("StokesCollage.jl", "w") do f
    write(f, string(symsim))
end

sim = include("StokesCollage.jl")
# sim = eval(symsim)
f = sim(sd, generate);

ω = eval_constant_primal_form(sd, SVector{3}([1.,0.,0.]));
u₀ = ComponentArray(v=ω,
                   P=ones(nv(sd))); # 17 
du = similar(u₀);
# constants_and_parameters=ComponentArray(μ=0.5,φ=zeros(length(ω)));
constants_and_parameters=ComponentArray(μ=1, lbP=lbP, rbP=rbP, tbP=tbP, bbP=bbP,φ=[1.0 for _ in 1:length(ω)])

# caps=ComponentArray(μ=1,φ=fill(1.0,length(ω)),boundaries=ComponentArray(lbP=lbP,rbP=rbP,tbP=tbP,bbP=bbP))

# f(du,u₀,constants_and_parameters,nothing) 
t=1e-2;
prob = ODEProblem(f,u₀,(0,t),constants_and_parameters);
soln = solve(prob, Tsit5());

function juxtapose(pode, msh, mshd, soln, vars; time::Float64=0)
    fig = Figure()
    foreach(enumerate(vars)) do (k,var)
        ax = CairoMakie.Axis(fig[1,2k-1], title=string(var))
        f = @match only(pode[incident(pode, var, :name), :type]) begin
            :Form1 => u -> norm.(♯(mshd, u, LLSDDSharp()))
            _ => identity
        end
        _m=mesh!(ax, msh, color=f(getproperty(soln(time), var)), colormap=:jet)
        Colorbar(fig[1,2k], _m)
    end
    fig
end
#
juxtapose(StokesCollage, s, sd, soln, [:P, :v], time=0.0) 

# fig = Figure()
# mesh(s, color=soln(t).P, colormap=:jet)
# fig = Figure()
# ax = CairoMakie.Axis(fig[1,1])
# msh = mesh!(ax, s, color=soln(t).P, colormap=:jet)
# Colorbar(fig[1,2], msh)
# fig

# mesh(s, color=soln(t).P, colormap=:jet)
# mesh(s, color=soln(t).v, colormap=:jet)
