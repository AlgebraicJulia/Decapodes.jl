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

using CUDA
using CUDA.CUSPARSE

## Original implementation can be found at: https://github.com/AlgebraicJulia/DECAPODES-Benchmarks

## Physics

begin 
  @info "Creating Physics"
  Diffusion = @decapode begin
    (T, Ṫ)::Form0
    ϕ::DualForm1
    # k::Constant

    # Fick's first law
    ϕ ==  ⋆₁(#=k=#0.4662225067399311 * d₀(T))
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
    # (kᵥ)::Constant
    V̇ == -(-(∧₁₀′(V, ⋆(d(V)))) + d(⋆(∧₁₁′(V, ⋆(V))))) + 
          #=kᵥ=# 1.716e-5 * (Δ₁(V) + (1/3) * (d₀(δ₁(V)))) / avg₀₁(ρ) +
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
    #(R₀)::Constant

    # ρ == p / (R₀ * T)
    ρ == p / (287.1015166027786 * T)
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

## Mesh Creation

s = EmbeddedDeltaSet2D("examples/diff_adv/su2_mesh_square_small_51.stl");
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Barycenter());
if ⋆(Val{1}, sd)[1,1] < 0.0
  orient_component!(s, 1, false)
end;


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

  wedge_mat(::Type{Val{(1,1)}}, s) = CuSparseMatrixCSC(2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s)))
  wedge_mat(::Type{Val{(2,0)}}, s) = CuSparseMatrixCSC(support_to_tri(s)*diag_vols(s)*tri_to_support(s))

  wedge_mat_11 = wedge_mat(Val{(1,1)},sd)
  wedge_mat_20 = wedge_mat(Val{(2,0)},sd)

  # wedge_edge(::Type{Val{(1,1)}}, s) = CuSparseMatrixCSC(dual_boundary(Val{2}, s) * 2.0 * (support_to_tri(s)*(diag_vols(s)*edge_to_support(s))))
  # wedge_edge(::Type{Val{(2,0)}}, s) = CuSparseMatrixCSC(dual_boundary(Val{2}, s) * support_to_tri(s)*(diag_vols(s)*tri_to_support(s)))

  dual_boundary_2 = CuSparseMatrixCSC{Float64}(dual_boundary(Val{2}, sd))
  wedge_edge(::Type{Val{(1,1)}}, s) = dual_boundary_2 * wedge_mat_11
  wedge_edge(::Type{Val{(2,0)}}, s) = dual_boundary_2 * wedge_mat_20

  function pd_wedge!(x, ::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat_11), caches=[CUDA.zeros(nv(s)), CUDA.zeros(ne(s))], kw...)
    broadcast!(*, caches[2], α, β)
    mul!(x, wedge_t[(1,1)], caches[2])
  end

  #= function pd_wedge(::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat_11), caches=[CUDA.zeros(nv(s)), CUDA.zeros(ne(s))], kw...)
    broadcast!(*, caches[2], α, β)
    wedge_t[(1,1)] * caches[2]
  end

  function pd_wedge(::Type{Val{(2,0)}}, s, α, β; wedge_t = Dict((2,0)=>wedge_mat_20), kw...)
    wedge_t[(2,0)] * (α .* β)
  end

  function pd_wedge(::Type{Val{(0,2)}}, s, α, β; wedge_t = nothing, kw...)
    α .* β
  end =#

  function init_wedge_ops(s)
    (d_mat=Dict(:d₁=>CuSparseMatrixCSC(d(Val{1}, s)), :dual_d₀=>CuSparseMatrixCSC(dual_derivative(Val{0}, s))),
    wedge_e=Dict((1,1)=>wedge_edge(Val{(1,1)},s), (2,0)=>wedge_edge(Val{(2,0)},s)),
    wedge_t=Dict((1,1)=>wedge_mat_11, (2,0)=>wedge_mat_20),
    caches=[CUDA.zeros(Float64, nv(s)), CUDA.zeros(Float64, ne(s)), CUDA.zeros(Float64, ntriangles(s))])
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
    CuSparseMatrixCSC(sparse(I,J,V, ne(s), ntriangles(s)*3))
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
    CuSparseMatrixCSC(sparse(I,J,V))
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
    CuSparseMatrixCSC(sparse(I,J,V))
  end

  function cp_2_1!(x, α, β, matrices)
    mul!(matrices[:α_cache], matrices[:t2c], α)
    mul!(matrices[:β_cache], matrices[:e2c], β)
    broadcast!(*, matrices[:β_cache], matrices[:α_cache], matrices[:β_cache])
    mul!(x, matrices[:cross], matrices[:β_cache])
    #x .= matrices[:cross] * ((matrices[:t2c]*α).*(matrices[:e2c]*β))
  end


  function avg_mat(::Type{Val{(0,1)}},s)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    for e in 1:ne(s)
        append!(J, [s[e,:∂v0],s[e,:∂v1]])
        append!(I, [e,e])
        append!(V, [0.5, 0.5])
    end
    CuSparseMatrixCSC(sparse(I,J,V))
  end
end

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

## BoundaryConditions

# Get boundaries of cylinders

begin
  locs = [(10.5, 0.0), (5.5, 0.0), (0.5, 0.0)]
  cyl = vcat(map(locs) do loc
      findall(p -> ( ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
  end...)
  fuzzy_bound = unique(vcat(incident(s, cyl, :∂v1)..., incident(s, cyl, :∂v0)...))
  cyl_edge = CuArray(filter(e -> (s[e, :∂v1] ∈ cyl)&&(s[e, :∂v0] ∈ cyl), fuzzy_bound))

  cyl_bound = vcat(map(locs) do loc
      findall(p -> (0.5^2 - 1e-4 <= ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
  end...)
  # fuzzy_boundary = unique(vcat(incident(s, cyl_bound, :∂v1)..., incident(s, cyl_bound, :∂v0)...))
  # cyl_bound_edge = filter(e -> (s[e, :∂v1] ∈ cyl_bound)&&(s[e, :∂v0] ∈ cyl_bound), fuzzy_boundary)
  cyl_inner = CuArray(filter(p -> !(p ∈ cyl_bound), cyl))
  # slip_edge = filter!(p -> !(p ∈ cyl_bound_edge), cyl_edge)

  k_col = CuArray(fill(k₁, ne(s)))
  k_col[cyl_edge] .= k₂
  # k = CuSparseMatrixCSC(diagm(k_col))
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
  # ∂₁ = boundary_inds(Val{1}, sd)

  ∂ₒ₀ = ∂₀[findall(p-> -9 <= p[1] <= 20 && -9 <= p[2] <= 9, s[∂₀, :point])]
  cu_∂ₒ₀ = CuArray(∂ₒ₀)
  lx = -10.0
  rx = 21.0
  ty = 15.5
  by = -15.5

  ∂ₗ₀ = ∂₀[findall(p-> p[1] <= lx + 1e-4, s[∂₀, :point])]
  ∂ᵣ₀ = ∂₀[findall(p-> p[1] >= rx - 1e-4, s[∂₀, :point])]
  ∂ₜ₀ = ∂₀[findall(p-> p[2] >= ty - 1e-4, s[∂₀, :point])]
  ∂ᵦ₀ = ∂₀[findall(p-> p[2] <= by + 1e-4, s[∂₀, :point])]
  ∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)

  # ∂ₗ₁ = bound_edges(s, ∂ₗ₀)
  # ∂ᵣ₁ = bound_edges(s, ∂ᵣ₀)
  # ∂ₑ₁ = bound_edges(s, ∂ₑ₀)

  # ∂₁₊ = adj_edges(s, ∂₀)

  # ∂ₗ₁₊ = adj_edges(s, ∂ₗ₀)
  # ∂ᵣ₁₊ = adj_edges(s, ∂ᵣ₀)
  ∂ₑ₁₊ = adj_edges(s, ∂ₑ₀)
  ∂_points = unique(vcat(s[∂ₑ₁₊, :∂v0], s[∂ₑ₁₊, :∂v1]))
  ∂ₑ₁₊ = bound_edges(s, ∂_points)
  cu_∂ₑ₁₊ = CuArray(∂ₑ₁₊)

  # ∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :∂v1], s[∂ₗ₁₊, :∂v0]))
  # ∂ᵣ₀₊ = unique(vcat(s[∂ᵣ₁₊, :∂v1], s[∂ᵣ₁₊, :∂v0]))
  ∂ₑ₀₊ = CuArray(unique(vcat(s[∂ₑ₁₊, :∂v1], s[∂ₑ₁₊, :∂v0])))

  c_objs = CuArray(fill(288.15, nv(s)))
  c_objs[cu_∂ₒ₀] .= 350.0
  velocity(p) = [3.402, 0.0, 0.0]
  gravity(p) = [0.0,0.0,0.0]
  v = CuArray(♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data);
  g = CuArray(♭(sd, DualVectorField(gravity.(sd[triangle_center(sd),:dual_point]))).data);
  p = CuArray([density for p in s[:point]] * (288.15 * R₀))
end
m_avg = avg_mat(Val{(0,1)}, sd)
# sim = eval(gensim(HeatXFer, code_target=cuda()))

wedge_cache = init_wedge_ops(sd)
v2comp = comp_support(sd);
cache_mat = Dict(:t2c => tri2comp(s, v2comp), :e2c => edge2comp(s, v2comp), :cross => changes(sd, v2comp),
                 :α_cache => CUDA.zeros(Float64, ntriangles(sd)*3), :β_cache => CUDA.zeros(Float64, ntriangles(sd)*3))
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :∂ₚ => begin 
    boundary = vcat(∂ₑ₀₊, cyl_inner)
      (x) -> begin
            # x[∂ₑ₀₊] .= 0
            # x[cyl_inner] .= 0
            x[boundary] .= 0
            x
      end
    end
    :avg₀₁ => x -> m_avg * x
    :∂ₜₐ => (x) -> begin
      x[cyl_inner] .= 0
      x
    end
    :∂ᵥ => begin 
    boundary = vcat(cyl_edge, cu_∂ₑ₁₊)
      (x) -> begin
        # x[cyl_edge] .= 0
        # x[cu_∂ₑ₁₊] .= 0
        x[boundary] .= 0
        x
      end
    end
    :∂ᵣ => begin 
    boundary = vcat(∂ₑ₀₊, cu_∂ₒ₀)
      (x) -> begin
        # x[∂ₑ₀₊] .= 0
        # x[cu_∂ₒ₀] .= 0
        x[boundary] .= 0
      end
    end
    #= :∧₁₀′ => (α, β) -> begin
      x = CUDA.zeros(Float64, ne(sd))
      cp_2_1!(x, β, α, cache_mat)
      x
    end
    =#
    :∧₁₀′ => (x, α, β) -> begin
      cp_2_1!(x, β, α, cache_mat)
    end
    #= :∧₁₁′ => (α, β) -> begin
      x = CUDA.zeros(Float64, nv(sd))
      pd_wedge!(x, Val{(1,1)}, sd, α, β; wedge_cache...)
      x
    =#
    :∧₁₁′ => (x, α, β) -> begin
      pd_wedge!(x, Val{(1,1)}, sd, α, β; wedge_cache...)
    end
    x => error("Unmatched operator $my_symbol")
  end
  return op
end

function simulate(mesh, operators, hodge = GeometricHodge())
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:577 =#
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:578 =#
  begin
      #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:174 =#
      (var"GenSim-M_d₁", d₁) = default_dec_cu_matrix_generate(mesh, :d₁, hodge)
      (var"GenSim-M_⋆₂", ⋆₂) = default_dec_cu_matrix_generate(mesh, :⋆₂, hodge)
      (var"GenSim-M_⋆₁", ⋆₁) = default_dec_cu_matrix_generate(mesh, :⋆₁, hodge)
      (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_cu_matrix_generate(mesh, :⋆₀⁻¹, hodge)
      (var"GenSim-M_d₀", d₀) = default_dec_cu_matrix_generate(mesh, :d₀, hodge)
      (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_cu_matrix_generate(mesh, :dual_d₁, hodge)
      (var"GenSim-M_⋆₀", ⋆₀) = default_dec_cu_matrix_generate(mesh, :⋆₀, hodge)
      (var"GenSim-M_dual_d₀", dual_d₀) = default_dec_cu_matrix_generate(mesh, :dual_d₀, hodge)     
      (var"GenSim-M_⋆₁⁻¹", ⋆₁⁻¹) = default_dec_cu_matrix_generate(mesh, :⋆₁⁻¹, hodge)
      (var"GenSim-M_∧₀₁", ∧₀₁) = default_dec_cu_matrix_generate(mesh, :∧₀₁, hodge)
      ∂ᵥ = operators(mesh, :∂ᵥ)
      ∂ᵣ = operators(mesh, :∂ᵣ)
      ∂ₚ = operators(mesh, :∂ₚ)
      # avg₀₁ = operators(mesh, :avg₀₁)
      avg₀₁ = avg_mat(Val{(0,1)}, sd)
      ∂ₜₐ = operators(mesh, :∂ₜₐ)
      (∧₁₀′) = operators(mesh, :∧₁₀′)
      (∧₁₁′) = operators(mesh, :∧₁₁′)
  end
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:579 =#
  begin
      #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:484 =#
      var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_d₀" * var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁"
      var"GenSim-ConMat_1" = (x->var"GenSim-M_GenSim-ConMat_1" * x)
      var"GenSim-M_GenSim-ConMat_2" = var"GenSim-M_dual_d₀" * var"GenSim-M_⋆₂" * var"GenSim-M_d₁"
      var"GenSim-ConMat_2" = (x->var"GenSim-M_GenSim-ConMat_2" * x)
      var"GenSim-M_GenSim-ConMat_3" = var"GenSim-M_⋆₂" * var"GenSim-M_d₁"
      var"GenSim-ConMat_3" = (x->var"GenSim-M_GenSim-ConMat_3" * x)
      var"GenSim-M_GenSim-ConMat_4" = var"GenSim-M_⋆₀⁻¹" * (var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁")          
      var"GenSim-ConMat_4" = (x->var"GenSim-M_GenSim-ConMat_4" * x)
      var"GenSim-M_GenSim-ConMat_5" = var"GenSim-M_d₀" * var"GenSim-M_⋆₀⁻¹"
      var"GenSim-ConMat_5" = (x->var"GenSim-M_GenSim-ConMat_5" * x)
      var"GenSim-M_GenSim-ConMat_6" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_⋆₀"
      var"GenSim-ConMat_6" = (x->var"GenSim-M_GenSim-ConMat_6" * x)
  end
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:580 =#
  begin
      #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:226 =#
      var"flow_•9" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•15" = CuVector{Float64}(undef, nparts(mesh, :E))
      Gensim_Var_5 = CuVector{Float64}(undef, nparts(mesh, :E))
      Gensim_Var_4 = CuVector{Float64}(undef, nparts(mesh, :E))
      Gensim_Var_1 = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•22" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•25" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•4" = CuVector{Float64}(undef, nparts(mesh, :Tri))
      var"diffusion_•3" = CuVector{Float64}(undef, nparts(mesh, :E))
      Gensim_Var_12 = CuVector{Float64}(undef, nparts(mesh, :V))
      Gensim_Var_15 = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•14" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•13" = CuVector{Float64}(undef, nparts(mesh, :E))
      Gensim_Var_13 = CuVector{Float64}(undef, nparts(mesh, :E))
      var"energy_•2" = CuVector{Float64}(undef, nparts(mesh, :V))
      ρ = CuVector{Float64}(undef, nparts(mesh, :V))
      Gensim_Var_16 = CuVector{Float64}(undef, nparts(mesh, :E))
      var"diffusion_•2" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•12" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•2" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•20" = CuVector{Float64}(undef, nparts(mesh, :V))
      Ṫ₁ = CuVector{Float64}(undef, nparts(mesh, :V))
      var"energy_•4" = CuVector{Float64}(undef, nparts(mesh, :V))
      flow_sum_3 = CuVector{Float64}(undef, nparts(mesh, :V))
      energy_Ṫₐ = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•6" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•19" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•24" = CuVector{Float64}(undef, nparts(mesh, :E))
      flow_sum_1 = CuVector{Float64}(undef, nparts(mesh, :E))
      flow_sum_2 = CuVector{Float64}(undef, nparts(mesh, :E))
      Ṫ = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•1" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•18" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•23" = CuVector{Float64}(undef, nparts(mesh, :E))
      ṗ = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•11" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•10" = CuVector{Float64}(undef, nparts(mesh, :E))
      V̇ = CuVector{Float64}(undef, nparts(mesh, :E))

      var"flow_•3" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•8" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•21" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"flow_•17" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"flow_•26" = CuVector{Float64}(undef, nparts(mesh, :E))

  end
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:581 =#
  f(du, u, p, t) = begin
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:581 =#
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:582 =#
          begin
              #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:251 =#
              V = u.V
              flow_G = u.flow_G
              T = u.T
              p = u.p
              flow_kᵥ = kᵥ
              energy_R₀ = R₀
              diffusion_k = k₂
              var"1" = 1.0
              var"3" = 3.0
              var"0.5" = 0.5
          end
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:583 =#
          mul!(var"flow_•9", var"GenSim-M_⋆₁", V)
          mul!(var"flow_•15", var"GenSim-M_GenSim-ConMat_1", V)
          mul!(Gensim_Var_5, var"GenSim-M_GenSim-ConMat_1", V)
          mul!(Gensim_Var_4, var"GenSim-M_GenSim-ConMat_2", V)
          var"GenSim-M_⋆₁⁻¹"(Gensim_Var_1, Gensim_Var_4)
          mul!(var"flow_•22", var"GenSim-M_⋆₁", V)
          mul!(var"flow_•25", var"GenSim-M_d₀", p)
          mul!(var"flow_•4", var"GenSim-M_GenSim-ConMat_3", V)
          mul!(var"diffusion_•3", var"GenSim-M_d₀", T)
          mul!(Gensim_Var_12, var"GenSim-M_GenSim-ConMat_6", p)
          mul!(Gensim_Var_15, var"GenSim-M_GenSim-ConMat_6", T)
          # var"flow_•3" = V ∧₁₀′ var"flow_•4"
          ∧₁₀′(var"flow_•3", V, var"flow_•4")
          # var"flow_•8" = V ∧₁₁′ var"flow_•9"
          ∧₁₁′(var"flow_•8", V, var"flow_•9")
          var"flow_•14" .= var"1" ./ var"3"
          var"flow_•13" .= var"flow_•14" .* var"flow_•15"
          # var"flow_•21" = V ∧₁₁′ var"flow_•22"
          ∧₁₁′(var"flow_•21", V, var"flow_•22")
          var"GenSim-M_∧₀₁"(Gensim_Var_13, Gensim_Var_12, V)
          var"energy_•2" .= energy_R₀ .* T
          ρ .= p ./ var"energy_•2"
          var"GenSim-M_∧₀₁"(Gensim_Var_16, Gensim_Var_15, V)
          var"diffusion_•2" .= diffusion_k .* var"diffusion_•3"
          var"flow_•12" .= (.+)(Gensim_Var_1, Gensim_Var_5)
          var"flow_•2" .= (.-)(var"flow_•3")
          # var"flow_•17" = avg₀₁(ρ)
          mul!(var"flow_•17", avg₀₁, ρ)
          mul!(var"flow_•20", var"GenSim-M_⋆₀⁻¹", var"flow_•21")
          # var"flow_•26" = avg₀₁(ρ)
          mul!(var"flow_•26", avg₀₁, ρ)
          mul!(Ṫ₁, var"GenSim-M_GenSim-ConMat_4", var"diffusion_•2")
          mul!(var"energy_•4", var"GenSim-M_GenSim-ConMat_4", Gensim_Var_16)
          mul!(flow_sum_3, var"GenSim-M_GenSim-ConMat_4", Gensim_Var_13)
          energy_Ṫₐ .= (.-)(var"energy_•4")
          energy_bc₀ = ∂ₜₐ(energy_Ṫₐ)
          mul!(var"flow_•6", var"GenSim-M_GenSim-ConMat_5", var"flow_•8")
          var"flow_•19" .= var"0.5" .* var"flow_•20"
          var"flow_•24" .= var"flow_•25" ./ var"flow_•26"
          flow_sum_1 .= (.+)(var"flow_•2", var"flow_•6")
          flow_sum_2 .= (.+)(var"flow_•12", var"flow_•13")
          Ṫ .= (.+)(energy_Ṫₐ, Ṫ₁)
          bcs_bc₀ = ∂ᵣ(Ṫ)
          var"flow_•1" .= (.-)(flow_sum_1)
          mul!(var"flow_•18", var"GenSim-M_d₀", var"flow_•19")
          var"flow_•23" .= (.-)(var"flow_•24")
          ṗ .= (.-)(flow_sum_3)
          var"flow_•11" .= flow_kᵥ .* flow_sum_2
          var"flow_•10" .= var"flow_•11" ./ var"flow_•17"
          V̇ .= (.+)(var"flow_•1", var"flow_•10", var"flow_•18", var"flow_•23", flow_G)
          bcs_bc₁ = ∂ᵥ(V̇)
          bcs_bc₀ = ∂ₚ(ṗ)
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:584 =#
          getproperty(du, :V) .= V̇
          getproperty(du, :p) .= ṗ
          getproperty(du, :T) .= Ṫ
      end
end
fₘ = simulate(sd, generate, GeometricHodge())

u₀ = ComponentArray(V=v, flow_G=g, T=c_objs, p=p)

constants_and_parameters = (energy_R₀=R₀, diffusion_k=k₂, flow_kᵥ=kᵥ)

tₑ = 2.0
tₑ = 0.5

dt = 0.01
@info "Solving ODE Problem"
prob = ODEProblem(fₘ, u₀, (0.0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5(), progress=true, progress_steps=1, dtmax=8e-4, saveat=dt, save_idxs=[:T, :p], p=g)

# fₘ(u₀, u₀, constants_and_parameters, 0)

densities(t) = Array(soln(t).p ./ (R₀ * soln(t).T))
begin
  frames = 150
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s, color=densities(0) , colormap=:jet, colorrange=extrema(densities(0)))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "CHT_cuda_3.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
      msh.color = densities(t)
  end
end
