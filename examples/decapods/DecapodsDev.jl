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

vect(s, e) = (s[s[e,:tgt], :point] - s[s[e,:src], :point]) * sign(1, s, e)
vect(s, e::AbstractVector) = [vect(s, el) for el in e]
t_vects(s,t) = vect(s, triangle_edges(s,t)) .* ([1,-1,1] * sign(2,s,t))
function comp_support(sd)
    vects = []
    for t in 1:ntriangles(sd)
        inds = triangle_edges(sd, t)
        for i in 1:3
            push!(vects, (t, inds[i]))
        end
    end
    v2comp = Dict{Tuple{Int64, Int64}, Int64}()
    for (i, v) in enumerate(vects)
        v2comp[v] = i
    end
    v2comp
end
function changes(s, v2comp)
    orient_vals = [1,-1,1]
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for t in 1:ntriangles(s)
        inds = triangle_edges(s, t)
        e_vects = t_vects(s,t)
        vals = zeros(1:3)
        for i in 1:3
            ns = [(i+1)%3 + 1, i%3+1]
            ort = e_vects[i] × (e_vects[i] × e_vects[ns[1]])
            n_ort = ort / norm(ort)
            append!(J, v2comp[(t,inds[i])])
            append!(I, inds[ns[1]])
            append!(V, dot(n_ort, e_vects[ns[1]]) * orient_vals[ns[1]] * sign(1, s, ns[1])* orient_vals[i]* sign(2,s,t) / 3.0)
            append!(J, v2comp[(t,inds[i])])
            append!(I, inds[ns[2]])
            append!(V, dot(n_ort, e_vects[ns[2]]) * orient_vals[ns[2]] * sign(1, s, ns[2])* orient_vals[i]* sign(2,s,t) / 3.0)
        end
    end
    sparse(I,J,V, ne(s), ntriangles(s)*3)
end
function edge2comp(s, v2comp)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for t in 1:ntriangles(sd)
        inds = triangle_edges(sd, t)
        for i in 1:3
            push!(I, v2comp[(t,inds[i])])
            push!(J, inds[i])
            push!(V, 1 / volume(Val{1}, s, inds[i]))
        end
    end
    sparse(I,J,V)
end
function tri2comp(s, v2comp)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for t in 1:ntriangles(sd)
        inds = triangle_edges(sd, t)
        for i in 1:3
            push!(I, v2comp[(t,inds[i])])
            push!(J, t)
            push!(V, 1)
        end
    end
    sparse(I,J,V)
end

function cp_2_1(α, β, matrices)
    matrices[:cross] * ((matrices[:t2c]*α).*(matrices[:e2c]*β))
    #matrices[:zeros]
end

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
  ∧₁₀′::Hom(Form1(X)⊗DualForm0(X), Form1(X))
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
  F1 ⋅ ∂ₜ(Form1(X)) == (F1 ⊗ (F2 ⋅ ⋆₂(X))) ⋅ ∧₁₀′ ⋅ neg₁
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


