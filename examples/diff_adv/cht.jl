using CombinatorialSpaces
using CombinatorialSpaces.DiscreteExteriorCalculus
using Catlab
using LinearAlgebra

using CairoMakie
using OrdinaryDiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Decapodes
using DiagrammaticEquations
using SparseArrays
using MLStyle
using ComponentArrays
using GeometryBasics
Point3D = Point3{Float64}

## Original implementation can be found at: https://github.com/AlgebraicJulia/DECAPODES-Benchmarks

## Physics

begin 
  @info "Creating Physics"
  Diffusion = @decapode begin
    (T, Ṫ)::Form0
    ϕ::DualForm1
    k::Constant

    # Fick's first law
    ϕ ==  ⋆₁(k * d₀(T))
    # Diffusion equation
    Ṫ == ⋆₀⁻¹(dual_d₁(ϕ))
  end

  Advection = @decapode begin
    (T, Ṫ)::Form0
    V::Form1
    Ṫ == -(⋆₀⁻¹(L(V, ⋆₀(T))))
  end

  Superposition = @decapode begin
    (Ṫ₁, Ṫ₂, Ṫ, T)::Form0
    Ṫ == Ṫ₁ + Ṫ₂
    ∂ₜ(T) == Ṫ
  end

  compose_diff_adv = @relation (T, Ṫ, V) begin
    diffusion(T, Ṫ₁)
    advection(T, Ṫ₂, V)
    superposition(Ṫ₁, Ṫ₂, Ṫ, T)
  end

  NavierStokes = @decapode begin
    (V, V̇, G)::Form1
    (T, ρ, ṗ, p)::Form0
    (kᵥ)::Constant
    V̇ == -(-(∧₁₀′(V, ⋆(d(V)))) + d(⋆(∧₁₁′(V, ⋆(V))))) + 
          kᵥ * (Δ₁(V) + (1/3) * (d₀(δ₁(V)))) / avg₀₁(ρ) +
          d₀(0.5 * ⋆(∧₁₁′(V, ⋆(V)))) +
          -(d₀(p) / avg₀₁(ρ)) +
          G
    ∂ₜ(V) == V̇
    ṗ == -(⋆₀⁻¹(L(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
    ∂ₜ(p) == ṗ
  end

  Energy = @decapode begin
    (V)::Form1
    (ρ, p, T, Ṫ, Ṫₐ, Ṫ₁, bc₀)::Form0
    (R₀)::Constant

    ρ == p / (R₀ * T)
    Ṫₐ == -(⋆₀⁻¹(L(V, ⋆₀(T))))
    ∂ₜₐ(Ṫₐ) == bc₀
    Ṫ == Ṫₐ + Ṫ₁
    ∂ₜ(T) == Ṫ
  end

  BoundaryConditions = @decapode begin 
    (V, V̇, bc₁)::Form1
    (Ṫ, ṗ, bc₀)::Form0 
    # no-slip edges
    ∂ᵥ(V̇) == bc₁
    # No change on left/right boundaries
    ∂ᵣ(Ṫ) == bc₀ ## Changed from partial t
    ∂ₚ(ṗ) == bc₀
  end

  compose_heat_xfer = @relation (V, ρ) begin
    flow(V, V̇, T, ρ, ṗ, p)
    energy(Ṫ, V, ρ, p, T, Ṫ₁)
    diffusion(T, Ṫ₁)
    bcs(Ṫ, ṗ, V, V̇)
  end

  @info "Composing Physics"
  HeatXfer_oapply = oapply(compose_heat_xfer,
                    [Open(NavierStokes, [:V, :V̇, :T, :ρ, :ṗ, :p]),
                    Open(Energy, [:Ṫ, :V, :ρ, :p, :T, :Ṫ₁]),
                    Open(Diffusion, [:T, :Ṫ]),
                    Open(BoundaryConditions, [:Ṫ, :ṗ, :V, :V̇])])

  @info "Extracting Conjugate Heat Transfer"
  HeatXFer = apex(HeatXfer_oapply)
  to_graphviz(HeatXFer)
end

## Special Ops

begin
  function edge_to_support(s)
    vals = Dict{Tuple{Int64, Int64}, Float64}()
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for e in 1:ne(s)
        de = elementary_duals(Val{1},s, e)
        ev = CombinatorialSpaces.volume(Val{1}, s, e)
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
        tv = CombinatorialSpaces.volume(Val{2}, s, t)
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

  diag_vols(s) = spdiagm([dual_volume(Val{2}, s, dt) for dt in 1:nparts(s, :DualTri)])

  wedge_mat(::Type{Val{(1,1)}}, s) = 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))
  wedge_mat(::Type{Val{(2,0)}}, s) = support_to_tri(s)*diag_vols(s)*tri_to_support(s)

  wedge_edge(::Type{Val{(1,1)}}, s) = dual_boundary(Val{2}, s) * 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))
  wedge_edge(::Type{Val{(2,0)}}, s) = dual_boundary(Val{2}, s) * support_to_tri(s)*diag_vols(s)*tri_to_support(s)

  function pd_wedge!(x, ::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), caches=[zeros(nv(s)), zeros(ne(s))], kw...)
    broadcast!(*, caches[2], α, β)
    mul!(x, wedge_t[(1,1)], caches[2])
  end

  function pd_wedge(::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), caches=[zeros(nv(s)), zeros(ne(s))], kw...)
    broadcast!(*, caches[2], α, β)
    wedge_t[(1,1)] * caches[2]
  end

  function pd_wedge(::Type{Val{(2,0)}}, s, α, β; wedge_t = Dict((2,0)=>wedge_mat(Val{(2,0)},s)), kw...)
    wedge_t[(2,0)] * (α .* β)
  end

  function pd_wedge(::Type{Val{(0,2)}}, s, α, β; wedge_t = nothing, kw...)
    α .* β
  end

  function init_wedge_ops(s)
    (d_mat=Dict(:d₁=>d(Val{1}, s), :dual_d₀=>dual_derivative(Val{0}, s)),
    wedge_e=Dict((1,1)=>wedge_edge(Val{(1,1)},s), (2,0)=>wedge_edge(Val{(2,0)},s)),
    wedge_t=Dict((1,1)=>wedge_mat(Val{(1,1)}, s), (2,0)=>wedge_mat(Val{(2,0)},s)),
    caches=[zeros(nv(s)), zeros(ne(s)), zeros(ntriangles(s))])
  end

  vect(s, e) = (s[s[e,:∂v1], :point] - s[s[e,:∂v0], :point]) * sign(1, s, e)
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
            push!(V, 1 / CombinatorialSpaces.volume(Val{1}, s, inds[i]))
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

  function cp_2_1!(x, α, β, matrices)
    mul!(matrices[:α_cache], matrices[:t2c], α)
    mul!(matrices[:β_cache], matrices[:e2c], β)
    broadcast!(*, matrices[:β_cache], matrices[:α_cache], matrices[:β_cache])
    mul!(x, matrices[:cross], matrices[:β_cache])
    #x .= matrices[:cross] * ((matrices[:t2c]*α).*(matrices[:e2c]*β))
  end
end

# TODO: Produce the implementations of these DEC operators

# # Lie derivative between two primal 1-forms
# Lie1Imp′ = @decapode ExtendedOperators begin
#   (F1, F1′, Ḟ1)::Form1{X}
#   Ḟ1 == i₀′(F1, d₁{X}(F1′)) + d₀{X}(i₁′(F1, F1′))
# end
# lie1_imp′ = diag2dwd(Lie1Imp′, in_vars = [:F1, :F1′], out_vars = [:Ḟ1])
# rules[:L₁′] = lie1_imp′

# # Internal product between a primal 1-form and a primal 2-form
# I0Imp′ = @decapode ExtendedOperators begin
#   (F1, Ḟ1)::Form1{X}
#   F2::Form2{X}
#   Ḟ1 == neg₁(∧₁₀′(F1, ⋆₂{X}(F2)))
# end
# i0_imp′ = diag2dwd(I0Imp′, in_vars = [:F1, :F2], out_vars = [:Ḟ1])
# rules[:i₀′] = i0_imp′

# # Internal product between two primal 1-forms
# I1Imp′ = @decapode ExtendedOperators begin
#   (F1, F1′)::Form1{X}
#   F0::Form0{X}
#   F0 == ⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
# end
# i1_imp′ = diag2dwd(I1Imp′, in_vars = [:F1, :F1′], out_vars = [:F0])
# rules[:i₁′] = i1_imp′

## Constants

begin
  cₚ = 1004.703 # Specific Heat at constant pressure
  kₜ = 0.0246295028571 #Thermal conductivity
  k_cyl = kₜ * 4
  kᵦ = 1.38064852e-23 # Boltzmann constant (m² kg/(s² K))

  density = 0.000210322
  R = kᵦ * 6.0221409e23 # kg⋅m²/(s²*K*mol)
  mol_mass = 28.96 # g/mol
  μ = 28.96 / 6.0221409e23 # mean molecular mass (g)
  R₀ = R / (mol_mass / 1000)

  # Heat diffusion constant in fluid
  k₁ = kₜ / (density * cₚ)

  # Heat diffusion constant in cylinder
  k₂ = k_cyl / (density * cₚ)

  kᵨ = 1e-3

  ν = 1.716e-5 # 0.0005081150545826

  kᵥ = ν# / density

  γ = 1 + R₀ / (cₚ - R₀) # Mayer's formula + definition of adiabatic constant
  e2t = 1/(cₚ * density)
  e2p = e2t * R₀
  t2e = 1/e2t
end
## Mesh Creation

s = EmbeddedDeltaSet2D("examples/diff_adv/su2_mesh_square_small_51.stl");
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Barycenter());
if ⋆(Val{1}, sd)[1,1] < 0.0
  orient_component!(s, 1, false)
end;

## BoundaryConditions

# Get boundaries of cylinders

begin
  locs = [(10.5, 0.0), (5.5, 0.0), (0.5, 0.0)]
  cyl = vcat(map(locs) do loc
      findall(p -> ( ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
  end...)
  fuzzy_bound = unique(vcat(incident(s, cyl, :∂v1)..., incident(s, cyl, :∂v0)...))
  cyl_edge = filter(e -> (s[e, :∂v1] ∈ cyl)&&(s[e, :∂v0] ∈ cyl), fuzzy_bound)

  cyl_bound = vcat(map(locs) do loc
      findall(p -> (0.5^2 - 1e-4 <= ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
  end...)
  fuzzy_boundary = unique(vcat(incident(s, cyl_bound, :∂v1)..., incident(s, cyl_bound, :∂v0)...))
  cyl_bound_edge = filter(e -> (s[e, :∂v1] ∈ cyl_bound)&&(s[e, :∂v0] ∈ cyl_bound), fuzzy_boundary)
  cyl_inner = filter(p -> !(p ∈ cyl_bound), cyl)
  slip_edge = filter!(p -> !(p ∈ cyl_bound_edge), cyl_edge)

  k_col = fill(k₁, ne(s))
  k_col[cyl_edge] .= k₂
  k = diagm(k_col)
end

# Get other boundaries
begin
  function boundary_inds(::Type{Val{1}}, s)
    collect(findall(x -> x != 0, boundary(Val{2},s) * fill(1,ntriangles(s))))
  end

  function boundary_inds(::Type{Val{0}}, s)
    ∂1_inds = boundary_inds(Val{1}, s)
    # TODO: Changed src and tgt to v0 and v1
    unique(vcat(s[∂1_inds,:∂v0],s[∂1_inds,:∂v1]))
  end

  function boundary_inds(::Type{Val{2}}, s)
    ∂1_inds = boundary_inds(Val{1}, s)
    inds = map([:∂e0, :∂e1, :∂e2]) do esym
      vcat(incident(s, ∂1_inds, esym)...)
    end
    unique(vcat(inds...))
  end

  function bound_edges(s, ∂₀)
    te = vcat(incident(s, ∂₀, :∂v1)...)
    se = vcat(incident(s, ∂₀, :∂v0)...)
    intersect(te, se)
  end

  function adj_edges(s, ∂₀)
    te = vcat(incident(s, ∂₀, :∂v1)...)
    se = vcat(incident(s, ∂₀, :∂v0)...)
    unique(vcat(te, se))
  end

  ∂₀ = boundary_inds(Val{0}, sd)
  ∂₁ = boundary_inds(Val{1}, sd)

  ∂ₒ₀ = ∂₀[findall(p-> -9 <= p[1] <= 20 && -9 <= p[2] <= 9, s[∂₀, :point])]

  lx = -10.0
  rx = 21.0
  ty = 15.5
  by = -15.5

  ∂ₗ₀ = ∂₀[findall(p-> p[1] <= lx + 1e-4, s[∂₀, :point])]
  ∂ᵣ₀ = ∂₀[findall(p-> p[1] >= rx - 1e-4, s[∂₀, :point])]
  ∂ₜ₀ = ∂₀[findall(p-> p[2] >= ty - 1e-4, s[∂₀, :point])]
  ∂ᵦ₀ = ∂₀[findall(p-> p[2] <= by + 1e-4, s[∂₀, :point])]
  ∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)

  ∂ₗ₁ = bound_edges(s, ∂ₗ₀)
  ∂ᵣ₁ = bound_edges(s, ∂ᵣ₀)
  ∂ₑ₁ = bound_edges(s, ∂ₑ₀)

  ∂₁₊ = adj_edges(s, ∂₀)

  ∂ₗ₁₊ = adj_edges(s, ∂ₗ₀)
  ∂ᵣ₁₊ = adj_edges(s, ∂ᵣ₀)
  ∂ₑ₁₊ = adj_edges(s, ∂ₑ₀)
  ∂_points = unique(vcat(s[∂ₑ₁₊, :∂v0], s[∂ₑ₁₊, :∂v1]))
  ∂ₑ₁₊ = bound_edges(s, ∂_points)

  ∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :∂v1], s[∂ₗ₁₊, :∂v0]))
  ∂ᵣ₀₊ = unique(vcat(s[∂ᵣ₁₊, :∂v1], s[∂ᵣ₁₊, :∂v0]))
  ∂ₑ₀₊ = unique(vcat(s[∂ₑ₁₊, :∂v1], s[∂ₑ₁₊, :∂v0]))

  c_objs = fill(288.15, nv(s))
  c_objs[∂ₒ₀] .= 350.0
  velocity(p) = [3.402, 0.0, 0.0]
  gravity(p) = [0.0,0.0,0.0]
  v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
  g = ♭(sd, DualVectorField(gravity.(sd[triangle_center(sd),:dual_point]))).data;
  p = [density for p in s[:point]] * (288.15 * R₀)
end

wedge_cache = init_wedge_ops(sd)
v2comp = comp_support(sd);
cache_mat = Dict(:t2c => tri2comp(s, v2comp), :e2c => edge2comp(s, v2comp), :cross => changes(sd, v2comp),
                 :α_cache => zeros(ntriangles(sd)*3), :β_cache => zeros(ntriangles(sd)*3))
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :∂ₚ => (x) -> begin
            x[∂ₑ₀₊] .= 0
            x[cyl_inner] .= 0
            x
          end
    :∂ₜₐ => (x) -> begin
      x[cyl_inner] .= 0
      x
    end
    :∂ᵥ => (x) -> begin
      x[cyl_edge] .= 0
      x[∂ₑ₁₊] .= 0
      x
    end
    :∂ᵣ => (x) -> begin
      x[∂ₑ₀₊] .= 0
      x[∂ₒ₀] .= 0
    end
    :∧₁₀′ => (α, β) -> begin
      x = zeros(ne(sd)) # TODO: Correct size?
      cp_2_1!(x, β, α, cache_mat)
      x
    end
    :∧₁₁′ => (α, β) -> begin
      x = zeros(nv(sd)) # TODO: Correct size?
      pd_wedge!(x, Val{(1,1)}, sd, α, β; wedge_cache...)
      x
    end
    x => error("Unmatched operator $my_symbol")
  end
  return op
end

fₘ = simulate(sd, generate)

u₀ = ComponentArray(V=v, flow_G=g, T=c_objs, p=p)

constants_and_parameters = (energy_R₀=R₀, diffusion_k=k₂, flow_kᵥ=kᵥ)

tₑ = 2.0
tₑ = 0.2

dt = 0.01
@info "Solving ODE Problem"
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5(), progress=true, progress_steps=1, dtmax=8e-4, saveat=dt, save_idxs=[:V, :flow_G, :T, :p], p=g)

# fₘ(u₀, u₀, constants_and_parameters, 0)

densities(t) = soln(t).p ./ (R₀ * soln(t).T)
inv_hdg_0 = inv_hodge_star(Val{0}, sd)
hdg_1 = ⋆(Val{1}, sd)
magnitudes(t) = sqrt.(abs.(1e-4 .+ inv_hdg_0*pd_wedge(Val{(1,1)}, sd, soln(t).V, hdg_1 * soln(t).V)))
begin
  frames = 90
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s, color=densities(0) , colormap=:seismic, colorrange=extrema(densities(0)))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "CHT_3.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
      msh.color = densities(t)
  end
end
