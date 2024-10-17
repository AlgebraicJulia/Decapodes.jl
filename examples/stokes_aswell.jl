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
StokesBoundaries = @decapode begin
    (lbP,rbP,tbP,bbP)::Constant
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
         :bbC => :bbP);
    restrictions=Dict(:lbP => :P, :rbP => :P, :tbP => :P, :bbP => :P))

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

# TODO: upstream restrict and mask to CombinatorialSpaces
"""restrict a form to a subset of the points.

- sd is the mesh,
- func is a function that chooses the submesh indices correspinding to the boundary
- form is the vector you want to restrict and 
"""
restrict(sd, func, form) = form[func(sd)]

restrict(indices, form) = form[indices]

"""mask a form to values on a subset of the points.

- sd is the mesh,
- form is the vector you want to restrict and 
- values is the vector you want to replace with
- func is a function that chooses the submesh indices correspinding to the boundary
"""
mask(sd::HasDeltaSet, func::Function, form, values) = begin
    #form[func(sd)] .= values
    restrict(sd, func, form) .= values
    form
end

mask(indices, form, values) = begin
    form[indices] .= values
    form
end

function generate(sd, my_symbol; hodge=DiagonalHodge())
    lwi = left_wall_idxs(sd)
    rwi = right_wall_idxs(sd)
    twi = top_wall_idxs(sd)
    bwi = bottom_wall_idxs(sd)
    op = @match my_symbol begin
        :restrict_lbP => ℓ -> restrict(lwi, ℓ)
        :restrict_rbP => ℓ -> restrict(rwi, ℓ)
        :restrict_tbP => ℓ -> restrict(twi, ℓ)
        :restrict_bbP => ℓ -> restrict(bwi, ℓ)
        :LeftBoundary => (ℓ,ℓ₀) -> mask(lwi, ℓ, ℓ₀)
        :RightBoundary => (ℓ,ℓ₀) -> mask(rwi, ℓ, ℓ₀)
        :TopBoundary => (ℓ,ℓ₀) -> mask(twi, ℓ, ℓ₀)
        :BottomBoundary => (ℓ,ℓ₀) -> mask(bwi, ℓ, ℓ₀)
        :.* => (x,y) -> x .* y
        _ => error("$my_symbol is not a supported operator.")
    end
    return (args...) -> op(args...)
end

lbP = fill(1.0, length(left_wall_idxs(sd)));
rbP = fill(1.0, length(right_wall_idxs(sd)));
tbP = fill(1.0, length(top_wall_idxs(sd)));
bbP = fill(1.0, length(bottom_wall_idxs(sd)));

# initialize data
symsim = gensim(StokesCollage);
open("StokesCollage.jl", "w") do f
    write(f, string(symsim))
end
#sim = include("StokesCollage.jl")
sim = eval(symsim)
f = sim(sd, generate);

ω = eval_constant_primal_form(sd, SVector{3}([1.,0.,0.]));
u₀ = ComponentArray(v=ω,
                   P=ones(nv(sd))); # 17 
du = similar(u₀);
constants_and_parameters=ComponentArray(μ=1,φ=[1.0 for _ in 1:length(ω)],lbP=lbP, rbP=rbP, tbP=tbP, bbP=bbP)

t=1e-5;
prob = ODEProblem(f,u₀,(0,t),constants_and_parameters);
soln = solve(prob, Tsit5());

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

juxtapose(StokesCollage, s, sd, soln, [:P, :v], time=t)
