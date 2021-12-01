# Our developed libraries
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus


using Catlab.Present
using Catlab.Graphics
using Catlab.CategoricalAlgebra
using LinearAlgebra

using Decapods.Simulations
using Decapods.Examples
using Decapods.Diagrams
using Decapods.Schedules

# Julia community libraries
using MeshIO
using CairoMakie
using Decapods.Debug
using DifferentialEquations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#####################################
# Primal-Primal wedge prototype def #
#####################################

using CombinatorialSpaces: volume
using SparseArrays

function edge_to_support(s)
    vals = Dict{Tuple{Int64, Int64}, Float64}()
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for e in 1:ne(s)
        de = elementary_duals(Val{1},s, e)
        ev = volume(Val{1}, s, e)
        dv = sum([dual_volume(Val{1}, s, d) for d in de])
        for d in de
            dt = incident(s, d, :D_∂e0)
            append!(I, dt)
            append!(J, fill(e, length(dt)))
            append!(V, fill(1/(dv*ev), length(dt)))
        end
    end
    sparse(I,J,V)
end

function tri_to_support(s)
    vals = Dict{Tuple{Int64, Int64}, Float64}()
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for t in 1:ntriangles(s)
        dt = vcat(incident(s, incident(s, triangle_center(s, t), :D_∂v0), :D_∂e1)...)
        tv = volume(Val{2}, s, t)
        append!(I, dt)
        append!(J, fill(t, length(dt)))
        append!(V, fill(1/tv#= * sign(Val{2}, s, t)=#, length(dt)))
    end
    sparse(I,J,V)
end

function support_to_tri(s)
    vals = Dict{Tuple{Int64, Int64}, Float64}()
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for t in 1:nv(s)
        dt = elementary_duals(Val{0},s, t)
        for d in dt
            push!(I, t)
            push!(J, d)
            push!(V, 1)
        end
    end
    sparse(I,J,V)
end

function support_edge_orient(s)
    # Dictionary of rel_o[edge][triangle] = sign
    rel_o = Vector{Dict{Int64, Int64}}(undef, ne(s))
    for o in 1:length(rel_o)
        rel_o[o] = Dict{Int64, Int64}()
    end
    for t in 1:ntriangles(s)
        edges = triangle_edges(s,t)
        rel_o[edges[1]][t] = sign(2,s,t) * sign(1,s,edges[1])
        rel_o[edges[2]][t] = #=-1 * =#sign(2,s,t) * sign(1,s,edges[2])
        rel_o[edges[3]][t] = sign(2,s,t) * sign(1,s,edges[3])
    end
    # Calculate edge for each support edge
    s2e = zeros(Int64, nparts(s, :DualTri))
    for e in 1:ne(s)
        de = elementary_duals(Val{1},s, e)
        for d in de
            dt = incident(s, d, :D_∂e0)
            s2e[dt] .= e
        end
    end
    # Calculate tri for each support edge
    s2t = zeros(Int64, nparts(s, :DualTri))
    for t in 1:ntriangles(s)
        dt = vcat(incident(s, incident(s, triangle_center(s, t), :D_∂v0), :D_∂e1)...)
        s2t[dt] .= t
    end
    # Calculate rel_o for each support edge
    s_orient = zeros(Int64, nparts(s, :DualTri))
    
    for i in 1:length(s2t)
        s_orient[i] = rel_o[s2e[i]][s2t[i]]
    end
    @show s_orient[1:10]
    s_orient
end

function support_to_edge(s)
    vals = Dict{Tuple{Int64, Int64}, Float64}()
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    s_orient = support_edge_orient(s)
    for e in 1:ne(s)
        de = elementary_duals(Val{1},s, e)
        for d in de
            dt = incident(s, d, :D_∂e0)
            append!(J, dt)
            append!(I, fill(e, length(dt)))
            append!(V, fill(1, length(dt)) .* s_orient[dt])
        end
    end
    sparse(I,J,V)
end

diag_vols(s) = spdiagm([dual_volume(Val{2}, s, dt) for dt in 1:nparts(s, :DualTri)])

wedge_mat(::Type{Val{(1,1)}}, s) = 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))
wedge_mat(::Type{Val{(2,0)}}, s) = support_to_tri(s)*diag_vols(s)*tri_to_support(s)

#wedge_edge(::Type{Val{(1,1)}}, s) = 0.5 * (support_to_edge(s)*diag_vols(s)*edge_to_support(s))
#wedge_edge(::Type{Val{(2,0)}}, s) = support_to_edge(s)*diag_vols(s)*tri_to_support(s)

wedge_edge(::Type{Val{(1,1)}}, s) = dual_boundary(Val{2}, s) * 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))
wedge_edge(::Type{Val{(2,0)}}, s) = dual_boundary(Val{2}, s) * support_to_tri(s)*diag_vols(s)*tri_to_support(s)

function pd_wedge(::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), kw...)
    wedge_t[(1,1)] * (α .* β)
end

function pd_wedge(::Type{Val{(2,0)}}, s, α, β; wedge_t = Dict((2,0)=>wedge_mat(Val{(2,0)},s)), kw...)
    wedge_t[(2,0)] * (α .* β)
end

function pd_wedge(::Type{Val{(0,2)}}, s, α, β; wedge_t = nothing, kw...)
    α .* β
end

# This might instead relate the calculated support volumes to the edges. It seems important that the resulting
# dimension is a dual 1
function pd_wedge(::Type{Val{(1,0)}}, s, α, β; d_mat = Dict(:d₁=>d(Val{1}, s), :dual_d₀=>dual_derivative(Val{0}, s)),
                                               wedge_e = Dict((1,1)=>wedge_edge(Val{(1,1)},s), 
                                                              (2,0)=>wedge_edge(Val{(2,0)},s)), kw...)
    #=(wedge_e[(2,0)] * ((d_mat[:d₁] * α) .* β)) -1 *=# (wedge_e[(1,1)] * (α .* (d_mat[:dual_d₀] * β)))
    zeros(ne(s))
end
function init_wedge_ops(s)
    (d_mat=Dict(:d₁=>d(Val{1}, s), :dual_d₀=>dual_derivative(Val{0}, s)),
     wedge_e=Dict((1,1)=>wedge_edge(Val{(1,1)},s), (2,0)=>wedge_edge(Val{(2,0)},s)),
     wedge_t=Dict((1,1)=>wedge_mat(Val{(1,1)}, s), (2,0)=>wedge_mat(Val{(2,0)},s)))
end



# Define used quantities
@present Flow2DQuantities(FreeExtCalc2D) begin
  X::Space
  C::Hom(munit(), Form0(X))     # concentration
  dC::Hom(munit(), Form0(X))     # concentration
  V::Hom(munit(), Form1(X))     # Flow field
  dV::Hom(munit(), Form1(X))    # Flow field time derivative
  P::Hom(munit(), DualForm1(X)) # Flow momentum
  p::Hom(munit(), Form0(X))     # pressure field
  dp::Hom(munit(), Form0(X))     # pressure field time derivative
  ϕ::Hom(munit(), DualForm1(X)) # negative diffusion flux
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  kᵥ::Hom(Form1(X), Form1(X))    # viscosity (usually scalar multiplication)
  kₚ::Hom(Form0(X), Form0(X))    # compressibility (usually scalar multiplication)
  m⁻¹::Hom(DualForm1(X), DualForm1(X))    # diffusivity (usually scalar multiplication)
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  ∂₀::Hom(Form0(X), Form0(X)) # all_wall boundary condition
  ∂₀₋::Hom(Form0(X), Form0(X)) # left/right wall boundary condition
  ∂₀₊::Hom(Form0(X), Form0(X)) # top/bottom wall boundary condition
  ∂₀ₛ::Hom(Form0(X), Form0(X)) # Sticking boundary condition
  ∂₁ₛ::Hom(Form1(X), Form1(X)) # Sticking boundary condition
  ∂₁ₗ₊::Hom(Form1(X), Form1(X)) # Sticking boundary condition
  ∂₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # In/Out edge flow boundary condition
  neg₁::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  dneg₀::Hom(DualForm0(X), DualForm0(X)) # In/Out edge flow boundary condition
  neg₀::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  mask₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  mask₀ₗ::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  half::Hom(Form0(X), Form0(X)) # half
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  bc₀::Hom(munit(), Form0(X))
  bc₁::Hom(munit(), Form1(X))
  bc₂::Hom(munit(), Form1(X))
  bc₃::Hom(munit(), Form0(X))
end

# Define Diffusion/Advection physics

@present DiffAdv2D <: Flow2DQuantities begin
  # Fick's first law
  ϕ == C ⋅ d₀(X) ⋅ k ⋅ ⋆₁(X)
  # Diffusion/advection equation
  dC == ϕ ⋅ dual_d₁(X) ⋅ ⋆₀⁻¹(X) + (V⊗(C⋅⋆₀(X))) ⋅ L₀ ⋅ ⋆₀⁻¹(X) ⋅ neg₀
  C ⋅ ∂ₜ(Form0(X)) == dC
  # Boundary condition
  C ⋅ ∂₀ == bc₀
end


@present IncompressibleFlow2D <: DiffAdv2D begin
  dV == (V ⊗ V) ⋅ L₁′ ⋅ neg₁ ⋅ mask₁ₑ + V ⋅ Δ₁(X) ⋅ kᵥ + (V ⊗ V) ⋅ i₁′ ⋅ half ⋅ d₀(X) + p ⋅ d₀(X) ⋅ neg₁
  V ⋅ ∂ₜ(Form1(X)) == dV
  dp == V ⋅ neg₁ ⋅ δ₁(X) ⋅ kₚ #+ (V ⊗ (p⋅⋆₀(X))) ⋅ L₀ ⋅ ⋆₀⁻¹(X) ⋅ mask₀ₗ ⋅ neg₀
  V ⋅ ∂₁ₛ == bc₁
  dC ⋅ ∂₀ₛ == bc₃
  dp ⋅ ∂₀₋ == bc₀
  dV ⋅ ∂₁ₗ₊ == bc₁
  p ⋅ ∂ₜ(Form0(X)) == dp
end;


## Special case operators for substitution
@present ExtendedOperators(FreeExtCalc2D) begin
  X::Space
  F0::Hom(munit(), Form0(X))
  F1::Hom(munit(), Form1(X))
  F2::Hom(munit(), Form2(X))
  dF0::Hom(munit(), DualForm0(X))
  dF1::Hom(munit(), DualForm1(X))
  dF2::Hom(munit(), DualForm2(X))
  neg::Hom(DualForm1(X), DualForm1(X)) # negative
  neg₁::Hom(Form1(X), Form1(X)) # negates 1-forms
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # negates 1-forms
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  i₀::Hom(Form1(X)⊗DualForm2(X), DualForm1(X))
  i₁::Hom(Form1(X)⊗DualForm1(X), DualForm0(X))
    
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  i₀′::Hom(Form1(X)⊗Form2(X), Form1(X))
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  ∧₁₁′::Hom(Form1(X)⊗DualForm1(X), DualForm2(X))
  ∧₁₀′::Hom(Form1(X)⊗DualForm0(X), DualForm1(X))
  F0′::Hom(munit(), Form0(X))
  F1′::Hom(munit(), Form1(X))
  F2′::Hom(munit(), Form2(X))
end

@present Rules <: ExtendedOperators begin
  i₀ == (id(Form1(X)) ⊗ ⋆₀⁻¹(X)) ⋅ ∧₁₀(X) ⋅ ⋆₁(X)
  i₁ == (id(Form1(X)) ⊗ ⋆₁⁻¹(X)) ⋅ ∧₁₁(X) ⋅ ⋆₂(X)
  L₀ == i₀ ⋅ dual_d₁(X)
  L₁ == (id(Form1(X))⊗dual_d₁(X)) ⋅ i₀ + i₁ ⋅ dual_d₀(X)
end;

@present Lie0Imp <: ExtendedOperators begin
  dF2 ⋅ ∂ₜ(DualForm2(X)) == (F1 ⊗ dF2) ⋅i₀ ⋅ dual_d₁(X)
end

@present Lie1Imp <: ExtendedOperators begin
  dF1 ⋅ ∂ₜ(DualForm1(X)) == (F1 ⊗ (dF1 ⋅ dual_d₁(X))) ⋅ i₀#= ⋅ dneg₁=# + (F1 ⊗ dF1) ⋅ i₁ ⋅ dual_d₀(X)
end

@present I0Imp <: ExtendedOperators begin
  dF1 ⋅ ∂ₜ(DualForm1(X)) == (F1 ⊗ (dF2 ⋅ ⋆₀⁻¹(X))) ⋅ ∧₁₀(X) ⋅ ⋆₁(X)
end

@present I1Imp <: ExtendedOperators begin
  dF0 ⋅ ∂ₜ(DualForm0(X)) == (F1 ⊗ (dF1 ⋅ ⋆₁⁻¹(X) ⋅ neg₁)) ⋅ ∧₁₁(X) ⋅ ⋆₂(X)
end

@present δ₁Imp <: ExtendedOperators begin
  F0 ⋅ ∂ₜ(Form0(X)) == F1 ⋅ ⋆₁(X) ⋅ dual_d₁(X) ⋅ ⋆₀⁻¹(X)
end

@present δ₂Imp <: ExtendedOperators begin
  F1 ⋅ ∂ₜ(Form1(X)) == F2 ⋅ ⋆₂(X) ⋅ dual_d₀(X) ⋅ ⋆₁⁻¹(X)
end

@present Δ0Imp <: ExtendedOperators begin
  F0 ⋅ ∂ₜ(Form0(X)) == F0 ⋅ d₀(X) ⋅ δ₁(X)
end

@present Δ1Imp <: ExtendedOperators begin
  F1 ⋅ ∂ₜ(Form1(X)) == F1 ⋅ (d₁(X) ⋅ δ₂(X) + δ₁(X) ⋅ d₀(X))
end

@present Lie1Imp′ <: ExtendedOperators begin
  F1 ⋅ ∂ₜ(Form1(X)) == (F1 ⊗ (F1′ ⋅ d₁(X))) ⋅ i₀′ + (F1 ⊗ F1′) ⋅ i₁′ ⋅ d₀(X)
end

@present I0Imp′ <: ExtendedOperators begin
  F1 ⋅ ∂ₜ(Form1(X)) == (F1 ⊗ (F2 ⋅ ⋆₂(X))) ⋅ ∧₁₀′ ⋅ ⋆₁⁻¹(X) ⋅ neg₁
end

@present I1Imp′ <: ExtendedOperators begin
  F0 ⋅ ∂ₜ(Form0(X)) == (F1 ⊗ (F1′ ⋅ ⋆₁(X))) ⋅ ∧₁₁′ ⋅ ⋆₀⁻¹(X)
end

lie0_imp_diag = eq_to_diagrams(Lie0Imp)
lie0_imp = diag2dwd(lie0_imp_diag)
tmp = lie0_imp.diagram[1, :outer_in_port_type]
lie0_imp.diagram[1, :outer_in_port_type] = lie0_imp.diagram[2, :outer_in_port_type]
lie0_imp.diagram[2, :outer_in_port_type] = tmp
tmp = lie0_imp.diagram[1, :in_src]
lie0_imp.diagram[1, :in_src] = lie0_imp.diagram[2, :in_src]
lie0_imp.diagram[2, :in_src] = tmp
to_graphviz(lie0_imp, orientation=LeftToRight)

lie1_imp_diag = eq_to_diagrams(Lie1Imp)
lie1_imp = diag2dwd(lie1_imp_diag)
tmp = lie1_imp.diagram[1, :outer_in_port_type]
lie1_imp.diagram[1, :outer_in_port_type] = lie1_imp.diagram[2, :outer_in_port_type]
lie1_imp.diagram[2, :outer_in_port_type] = tmp
lie1_imp.diagram[1, :in_src] = 2
lie1_imp.diagram[2, :in_src] = 1
lie1_imp.diagram[3, :in_src] = 1
lie1_imp.diagram[4, :in_src] = 2
to_graphviz(lie1_imp, orientation=LeftToRight)

i0_imp_diag = eq_to_diagrams(I0Imp)
i0_imp = diag2dwd(i0_imp_diag)
rem_part!(i0_imp.diagram, :OuterInPort, 1)

#tmp = i0_imp.diagram[1, :in_src]
#i0_imp.diagram[1, :in_src] = i0_imp.diagram[2, :in_src]
#i0_imp.diagram[2, :in_src] = tmp

to_graphviz(i0_imp, orientation=LeftToRight)

i1_imp_diag = eq_to_diagrams(I1Imp)
i1_imp = diag2dwd(i1_imp_diag)
rem_part!(i1_imp.diagram, :OuterInPort, 1)

#tmp = i1_imp.diagram[1, :in_src]
#i1_imp.diagram[1, :in_src] = i1_imp.diagram[2, :in_src]
#i1_imp.diagram[2, :in_src] = tmp

to_graphviz(i1_imp, orientation=LeftToRight)

δ₁_imp_diag = eq_to_diagrams(δ₁Imp)
δ₁_imp = diag2dwd(δ₁_imp_diag)
rem_part!(δ₁_imp.diagram, :OuterInPort, 1)

to_graphviz(δ₁_imp, orientation=LeftToRight)

δ₂_imp_diag = eq_to_diagrams(δ₂Imp)
δ₂_imp = diag2dwd(δ₂_imp_diag)
rem_part!(δ₂_imp.diagram, :OuterInPort, 1)

to_graphviz(δ₂_imp, orientation=LeftToRight)

Δ0_imp_diag = eq_to_diagrams(Δ0Imp)
Δ0_imp = diag2dwd(Δ0_imp_diag)

to_graphviz(Δ0_imp, orientation=LeftToRight)

Δ1_imp_diag = eq_to_diagrams(Δ1Imp)
Δ1_imp = diag2dwd(Δ1_imp_diag)

to_graphviz(Δ1_imp, orientation=LeftToRight)

lie1_imp_diag′ = eq_to_diagrams(Lie1Imp′)
lie1_imp′ = diag2dwd(lie1_imp_diag′)
#=tmp = lie1_imp′.diagram[1, :outer_in_port_type]
lie1_imp.diagram[1, :outer_in_port_type] = lie1_imp.diagram[2, :outer_in_port_type]
lie1_imp.diagram[2, :outer_in_port_type] = tmp
lie1_imp.diagram[1, :in_src] = 2
lie1_imp.diagram[2, :in_src] = 1
lie1_imp.diagram[3, :in_src] = 1
lie1_imp.diagram[4, :in_src] = 2=#
to_graphviz(lie1_imp′, orientation=LeftToRight)

i0_imp_diag′ = eq_to_diagrams(I0Imp′)
i0_imp′ = diag2dwd(i0_imp_diag′)
#rem_part!(i0_imp.diagram, :OuterInPort, 1)

#tmp = i0_imp.diagram[1, :in_src]
#i0_imp.diagram[1, :in_src] = i0_imp.diagram[2, :in_src]
#i0_imp.diagram[2, :in_src] = tmp

to_graphviz(i0_imp′, orientation=LeftToRight)

i1_imp_diag′ = eq_to_diagrams(I1Imp′)
i1_imp′ = diag2dwd(i1_imp_diag′)
rem_part!(i1_imp′.diagram, :OuterInPort, 1)

#tmp = i1_imp.diagram[1, :in_src]
#i1_imp.diagram[1, :in_src] = i1_imp.diagram[2, :in_src]
#i1_imp.diagram[2, :in_src] = tmp

to_graphviz(i1_imp′, orientation=LeftToRight)

p_diag = eq_to_diagrams(IncompressibleFlow2D)
to_graphviz(p_diag ; edge_len="1.3")
#Graph(diag)

(new_dwd =diag2dwd(p_diag, calc_states=[], clean=true))
to_graphviz(new_dwd, orientation=LeftToRight)

(exp_dwd = Examples.expand_dwd(new_dwd, Dict(:L₀ => lie0_imp, :i₀ => i0_imp, :L₁ => lie1_imp, :i₁ => i1_imp,
                                             :δ₁ => δ₁_imp, :δ₂ => δ₂_imp, :Δ₀ => Δ0_imp, :Δ₁ => Δ1_imp,
                                             :i₀′ => i0_imp′, :L₁′ => lie1_imp′, :i₁′ => i1_imp′)))
to_graphviz(exp_dwd, orientation=LeftToRight)

s = EmbeddedDeltaSet2D("../../../../../meshes/cylinder_forward_inlet_outlet_2.stl")
orient_component!(s, 1, false)
sd = dual(s, subdiv=Barycenter());

k = 0.1
kₚ = 1600.0
kᵥ = 2.0
pressure = -1.0# -5
ramp_up = 10.0#0.1
funcs = sym2func(sd)
∂₀ = Examples.boundary_inds(Val{0}, s)
∂₁ = Examples.boundary_inds(Val{1}, s)

∂ₗ₀ = ∂₀[findall(p-> p[1] <= -50.0, s[∂₀, :point])]
∂ᵣ₀ = ∂₀[findall(p-> p[1] >= 50.0, s[∂₀, :point])]
∂ₜ₀ = ∂₀[findall(p-> p[2] >= 15.0, s[∂₀, :point])]
∂ᵦ₀ = ∂₀[findall(p-> p[2] <= -15.0, s[∂₀, :point])]
∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)
∂ₒ₀ = ∂₀[findall(p-> -49 <= p[1] <= 49 && -14 <= p[2] <= 14, s[∂₀, :point])]

∂ₗ₁ = Examples.bound_edges(s, ∂ₗ₀)
∂ᵣ₁ = Examples.bound_edges(s, ∂ᵣ₀)
∂ₗ₁₊ = Examples.adj_edges(s, ∂ₗ₀)
∂ᵣ₁₊ = Examples.adj_edges(s, ∂ᵣ₀)
∂ₒ₁ = Examples.bound_edges(s, ∂ₒ₀)


∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :tgt], s[∂ₗ₁₊, :src]))

funcs[:mcopy] = Dict(:operator => I(ne(sd)), :type => MatrixFunc())
funcs[:k] = Dict(:operator => k * I(ne(sd)), :type => MatrixFunc())
funcs[:half] = Dict(:operator => 0.5 * I(nv(sd)), :type => MatrixFunc())
funcs[:kₚ] = Dict(:operator => kₚ * I(nv(sd)), :type => MatrixFunc())
funcs[:kᵥ] = Dict(:operator => kᵥ * I(ne(sd)), :type => MatrixFunc())
funcs[:dneg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₀] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())
funcs[:∂₀] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀] .= 0), :type => InPlaceFunc())

ramp_up_pressure(x′, x, t) = begin
    x′ .= x
    if t < (-1 * pressure * 50) / ramp_up
        x′[∂ₗ₀] .= ramp_up
        x′[∂ᵣ₀] .= -(ramp_up)
    else
        x′[∂ₗ₀] .= 0
        x′[∂ᵣ₀] .= 0
    end
end

funcs[:∂₀₋] = Dict(:operator => ramp_up_pressure, :type => TDInPlaceFunc())
funcs[:∂ₗ₀₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₗ₊] = Dict(:operator => (x′,x) -> (x′ .= x), :type => InPlaceFunc())
#funcs[:∂₁ₗ₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0), :type => InPlaceFunc())
funcs[:∂₀ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₀] .= 0; x′[∂ₑ₀] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₁] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₀ₗ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())

wedge_cache = init_wedge_ops(sd)
funcs[:∧₁₀′] = Dict(:operator => (x′, α, β) -> (x′ .= pd_wedge(Val{(1,0)},s, α, β; wedge_cache...)), :type => InPlaceFunc())
funcs[:∧₁₁′] = Dict(:operator => (x′, α, β) -> (x′ .= pd_wedge(Val{(1,1)},s, α, β; wedge_cache...)), :type => InPlaceFunc())

const hodge_LU = lu(⋆(Val{1}, sd))
funcs[:⋆₁⁻¹] = Dict(:operator => (x′,x) -> (x′ .= -1 * (hodge_LU \ x)), :type => InPlaceFunc())

c_objs = zeros(nv(s))
c_objs[∂ₒ₀] .= 1e-3
#c_objs = [sum(pdf(c_dist, [p[1] - cx, p[2]]) for cx in [0.5, 5.5, 10.5] ) for p in s[:point]]
#velocity(p) = [30.0 * exp(-(p[1]+50)/1),0.0,0.0]
velocity(p) = [0.0,0.0,0.0]
#velocity(p) = [10 * (1-(p[2]/15)^2) * exp(-(p[1]+50)/5),0.0,0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
p = [(-1 * (p[1]-50))/100.0 * 30.0 for p in s[:point]]
#p = [0.0 for p in s[:point]]
u0 = vcat(c_objs, v, p)
GC.gc()
@show ⋆(Val{1},sd)[1,1]
fig, ax, ob = wireframe(s)
ax.aspect = AxisAspect(3)
fig

#=s = EmbeddedDeltaSet2D("../../meshes/su2_mesh_square_coarse.stl")
orient_component!(s, 1, false)
sd = dual(s);

k = 0.1
kₚ = 1600.0
kᵥ = 2.0
pressure = -1.0# -5
ramp_up = 10.0#0.1
funcs = sym2func(sd)
∂₀ = Examples.boundary_inds(Val{0}, s)
∂₁ = Examples.boundary_inds(Val{1}, s)

∂ₗ₀ = ∂₀[findall(p-> p[1] <= -30.0, s[∂₀, :point])]
∂ᵣ₀ = ∂₀[findall(p-> p[1] >= 41.0, s[∂₀, :point])]
∂ₜ₀ = ∂₀[findall(p-> p[2] >= 35.5, s[∂₀, :point])]
∂ᵦ₀ = ∂₀[findall(p-> p[2] <= -35.5, s[∂₀, :point])]
∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)
∂ₒ₀ = ∂₀[findall(p-> -20 <= p[1] <= 20 && -20 <= p[2] <= 20, s[∂₀, :point])]

∂ₗ₁ = Examples.bound_edges(s, ∂ₗ₀)
∂ᵣ₁ = Examples.bound_edges(s, ∂ᵣ₀)
∂ₗ₁₊ = Examples.adj_edges(s, ∂ₗ₀)
∂ᵣ₁₊ = Examples.adj_edges(s, ∂ᵣ₀)
∂ₒ₁ = Examples.bound_edges(s, ∂ₒ₀)


∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :tgt], s[∂ₗ₁₊, :src]))

funcs[:mcopy] = Dict(:operator => I(ne(sd)), :type => MatrixFunc())
funcs[:k] = Dict(:operator => k * I(ne(sd)), :type => MatrixFunc())
funcs[:half] = Dict(:operator => 0.5 * I(nv(sd)), :type => MatrixFunc())
funcs[:kₚ] = Dict(:operator => kₚ * I(nv(sd)), :type => MatrixFunc())
funcs[:kᵥ] = Dict(:operator => kᵥ * I(ne(sd)), :type => MatrixFunc())
funcs[:dneg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₀] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())
funcs[:∂₀] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀] .= 0), :type => InPlaceFunc())

ramp_up_pressure(x′, x, t) = begin
    x′ .= x
    if t < (-1 * pressure * 50) / ramp_up
        x′[∂ₗ₀] .= ramp_up
        x′[∂ᵣ₀] .= -(ramp_up)
    else
        x′[∂ₗ₀] .= 0
        x′[∂ᵣ₀] .= 0
    end
end

funcs[:∂₀₋] = Dict(:operator => ramp_up_pressure, :type => TDInPlaceFunc())
funcs[:∂ₗ₀₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₗ₊] = Dict(:operator => (x′,x) -> (x′ .= x), :type => InPlaceFunc())
#funcs[:∂₁ₗ₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0), :type => InPlaceFunc())
funcs[:∂₀ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₀] .= 0; x′[∂ₑ₀] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₁] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₀ₗ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())

wedge_cache = init_wedge_ops(sd)
funcs[:∧₁₀′] = Dict(:operator => (x′, α, β) -> (x′ .= pd_wedge(Val{(1,0)},s, α, β; wedge_cache...)), :type => InPlaceFunc())
funcs[:∧₁₁′] = Dict(:operator => (x′, α, β) -> (x′ .= pd_wedge(Val{(1,1)},s, α, β; wedge_cache...)), :type => InPlaceFunc())

const hodge_LU = lu(⋆(Val{1}, sd))
funcs[:⋆₁⁻¹] = Dict(:operator => (x′,x) -> (x′ .= -1 * (hodge_LU \ x)), :type => InPlaceFunc())

#Gaussian curve
c_objs = zeros(nv(s))
c_objs[∂ₒ₀] .= 1e-3
#c_objs = [sum(pdf(c_dist, [p[1] - cx, p[2]]) for cx in [0.5, 5.5, 10.5] ) for p in s[:point]]
velocity(p) = [3.14 * exp(-(p[1] + 30)/1),0.0,0.0]
velocity(p) = [0.0,0.0,0.0]
#velocity(p) = [10 * (1-(p[2]/15)^2) * exp(-(p[1]+50)/5),0.0,0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
p = [(-1 * (p[1]-50))/100.0 * 30.0 for p in s[:point]]
#p = [0.0 for p in s[:point]]
u0 = vcat(c_objs, v, p)
GC.gc()
@show ⋆(Val{1},sd)[1,1]
fig, ax, ob = wireframe(s)
ax.aspect = AxisAspect(1.0)
fig=#


new_dwd = deepcopy(exp_dwd)
Examples.contract_matrices!(new_dwd, funcs)

to_graphviz(new_dwd, orientation=LeftToRight)


func, _ = gen_sim(new_dwd, funcs, sd; autodiff=false);

prob = ODEProblem(func, u0, (0, 10.0))
sol1 = solve(prob, Tsit5(), progress=true, progress_steps=100, p=v, saveat=0.1)

times = range(0, sol1.t[end], length=150)
colors = [sol1(t)#=[(end-nv(s)+1):end]=# for t in times]

figure, axis, scatter_thing = mesh(s, color=colors[1],
                                   colorrange=(0,1e-3))#minimum(vcat(colors...)),maximum(vcat(colors...))))
axis.aspect = AxisAspect(100.0/30.0)
framerate = 30

record(figure, "conc_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
  scatter_thing.color = colors[i]
end

times = range(0, sol1.t[end], length=150)
colors = [sol1(t)[(end-nv(s)+1):end] for t in times]

figure, axis, scatter_thing = mesh(s, color=colors[1],
                                   colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
axis.aspect = AxisAspect(100.0/30.0)
framerate = 30

record(figure, "press_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
  scatter_thing.color = colors[i]
end
