# Our developed libraries
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using CombinatorialSpaces: volume


using Catlab.Present
using Catlab.Graphics
using Catlab.CategoricalAlgebra

using Decapods.Simulations
using Decapods.Examples
using Decapods.Diagrams
using Decapods.Schedules

# Julia community libraries
using LinearAlgebra
using MeshIO
using CairoMakie
using DifferentialEquations
using Logging: global_logger
using SparseArrays
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Decapods.Debug
include("DecapodsDev.jl")

#####################################
# Physical Quantities and Equations #
#####################################

@present Flow2DQuantities(FreeExtCalc2D) begin
  X::Space
  
  # Physical variables
  T::Hom(munit(), Form0(X))     # temp
  dT::Hom(munit(), Form0(X))     # change in temp
  ρ::Hom(munit(), Form0(X))     # density
  dρ::Hom(munit(), Form0(X))     # change in density
  V::Hom(munit(), Form1(X))     # Flow field
  dV::Hom(munit(), Form1(X))    # Flow field time derivative
  P::Hom(munit(), DualForm1(X)) # Flow momentum
  p::Hom(munit(), Form0(X))     # pressure field
  dp::Hom(munit(), Form0(X))     # pressure field time derivative
  ϕ::Hom(munit(), DualForm1(X)) # negative diffusion flux
  
  # Multiplicative constants
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  R₀::Hom(Form0(X), Form0(X))    # Ideal gas constant (usually scalar multiplication)
  kᵥ::Hom(Form1(X), Form1(X))    # viscosity (usually scalar multiplication)
  kᵥₚ::Hom(Form0(X), Form0(X))    # viscosity (usually scalar multiplication)
  kₚ::Hom(Form0(X), Form0(X))    # compressibility (usually scalar multiplication)
 
  # Boundary conditions
  ∂₀::Hom(Form0(X), Form0(X))   # all_wall boundary condition
  ∂₀₋::Hom(Form0(X), Form0(X))  # left/right wall boundary condition
  ∂₀₊::Hom(Form0(X), Form0(X))  # top/bottom wall boundary condition
  ∂₀ₛ::Hom(Form0(X), Form0(X))  # Sticking boundary condition
  ∂₁ₛ::Hom(Form1(X), Form1(X))  # Sticking boundary condition
  ∂₁ₗ₊::Hom(Form1(X), Form1(X)) # Left Edge boundary condition
  ∂₁ₑ::Hom(Form1(X), Form1(X))    # In/Out edge flow boundary condition
  mask₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  mask₀ₗ::Hom(Form0(X), Form0(X)) # Left edge boundary condition
  mask₀ₛ::Hom(Form0(X), Form0(X)) # Flow restriction within cylinder
  bc₀::Hom(munit(), Form0(X))
  bc₁::Hom(munit(), Form1(X))
  bc₂::Hom(munit(), Form1(X))
  bc₃::Hom(munit(), Form0(X)) 

  # Operators
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  half::Hom(Form0(X), Form0(X))     # multiplication by 1/2
  third::Hom(Form1(X), Form1(X))    # multiplication by 1/3
  mult₀::Hom(Form0(X)⊗Form0(X), Form0(X))# elementwise multiplication
  avg₀₁::Hom(Form0(X), Form1(X))         # Interpolation from 0-forms to 1-forms
  div₁::Hom(Form1(X)⊗Form1(X), Form1(X)) # elementwise division
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # Negation of dual 1-form 
  neg₁::Hom(Form1(X), Form1(X))          # Negation of 1-form 
  dneg₀::Hom(DualForm0(X), DualForm0(X)) # Negation of dual 0-forms 
  neg₀::Hom(Form0(X), Form0(X))          # Negation of 0-forms 
end

@present DiffAdv2D <: Flow2DQuantities begin
  # Fick's first law
  ϕ == T ⋅ d₀(X) ⋅ k ⋅ ⋆₁(X)
  # Diffusion/advection equation
  dT == ϕ ⋅ dual_d₁(X) ⋅ ⋆₀⁻¹(X) + (V⊗(T⋅⋆₀(X))) ⋅ L₀ ⋅ ⋆₀⁻¹(X) ⋅ neg₀ ⋅ mask₀ₛ ⋅ mask₀ₗ
  T ⋅ ∂ₜ(Form0(X)) == dT
  # Boundary condition
  T ⋅ ∂₀ == bc₀
end
@present OneWayCoupledFlow2D <: DiffAdv2D begin
  dV == (V ⊗ V) ⋅ L₁′ ⋅ neg₁ ⋅ mask₁ₑ + V ⋅ (Δ₁(X) + δ₁(X)⋅d₀(X)⋅third) ⋅ kᵥ + (V ⊗ V) ⋅ i₁′ ⋅ half ⋅ d₀(X) + p ⋅ d₀(X) ⋅ neg₁
  V ⋅ ∂ₜ(Form1(X)) == dV
  dp == (V ⋅ neg₁ ⋅ δ₁(X) ⋅ kₚ + p ⋅ Δ₀(X) ⋅ kᵥₚ) ⋅ mask₀ₛ
  V ⋅ ∂₁ₛ == bc₁
  dT ⋅ ∂₀ₛ == bc₃
  dp ⋅ ∂₀₋ == bc₀
  dV ⋅ ∂₁ₗ₊ == bc₁
  dV ⋅ ∂₁ₛ == bc₁
  p ⋅ ∂ₜ(Form0(X)) == dp
end;

@present TwoWayCoupledFlow2D <: DiffAdv2D begin
  dV == (V ⊗ V) ⋅ L₁′ ⋅ neg₁ ⋅ mask₁ₑ + V ⋅ (Δ₁(X) + δ₁(X)⋅d₀(X)⋅third) ⋅ kᵥ + (V ⊗ V) ⋅ i₁′ ⋅ half ⋅ d₀(X) + ((p ⋅ d₀(X))⊗(ρ ⋅ avg₀₁)) ⋅ div₁ ⋅ neg₁
  V ⋅ ∂ₜ(Form1(X)) == dV
  ρ ⋅ ∂ₜ(Form0(X)) == (V ⊗ (ρ⋅⋆₀(X))) ⋅ L₀ ⋅ ⋆₀⁻¹(X) ⋅ mask₀ₛ ⋅ mask₀ₗ ⋅ neg₀
  p == (T ⊗ ρ) ⋅ mult₀ ⋅ R₀
  V ⋅ ∂₁ₛ == bc₁
  dT ⋅ ∂₀ₛ == bc₃
  dV ⋅ ∂₁ₗ₊ == bc₁
  dV ⋅ ∂₁ₛ == bc₁
end;

# Generate equational diagram
p_diag = eq_to_diagrams(OneWayCoupledFlow2D)

# Schedule equational diagram
new_dwd =diag2dwd(p_diag, calc_states=[], clean=true)

# Provide complex operator definitions
exp_dwd = Examples.expand_dwd(new_dwd, Dict(:L₀ => lie0_imp, :i₀ => i0_imp, :L₁ => lie1_imp, :i₁ => i1_imp,
                                             :δ₁ => δ₁_imp, :δ₂ => δ₂_imp, :Δ₀ => Δ0_imp, :Δ₁ => Δ1_imp,
                                             :i₀′ => i0_imp′, :L₁′ => lie1_imp′, :i₁′ => i1_imp′))


#########################################
# Import Mesh and Define BCs/constants #
#########################################
s1 = EmbeddedDeltaSet2D("../../../../meshes/su2_mesh_square_small_51.stl")
s = EmbeddedDeltaSet2D{Bool, Point{3, Float64}}()
copy_parts!(s, s1)
sd = dual(s);

# Extract boundaries of cylinders
locs = [(10.5, 0.0), (5.5, 0.0), (0.5, 0.0)]
cyl = vcat(map(locs) do loc
    findall(p -> ( ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
end...)
fuzzy_bound = unique(vcat(incident(s, cyl, :tgt)..., incident(s, cyl, :src)...))
cyl_edge = filter(e -> (s[e, :tgt] ∈ cyl)&&(s[e, :src] ∈ cyl), fuzzy_bound)
cyl_bound = vcat(map(locs) do loc
    findall(p -> (0.5^2 - 1e-4 <= ((p[1] - loc[1])^2 + (p[2] - loc[2])^2) <= 0.5^2 + 1e-4),s[:point])
end...)
fuzzy_boundary = unique(vcat(incident(s, cyl_bound, :tgt)..., incident(s, cyl_bound, :src)...))
cyl_bound_edge = filter(e -> (s[e, :tgt] ∈ cyl_bound)&&(s[e, :src] ∈ cyl_bound), fuzzy_boundary)
cyl_inner = filter(p -> !(p ∈ cyl_bound), cyl)

∂₀ = Examples.boundary_inds(Val{0}, s)
∂₁ = Examples.boundary_inds(Val{1}, s)

∂ₗ₀ = ∂₀[findall(p-> p[1] <= -10.0, s[∂₀, :point])]
∂ᵣ₀ = ∂₀[findall(p-> p[1] >= 21.0, s[∂₀, :point])]
∂ₜ₀ = ∂₀[findall(p-> p[2] >= 15.5, s[∂₀, :point])]
∂ᵦ₀ = ∂₀[findall(p-> p[2] <= -15.5, s[∂₀, :point])]
∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)
∂ₒ₀ = ∂₀[findall(p-> -9 <= p[1] <= 20 && -9 <= p[2] <= 9, s[∂₀, :point])]

∂ₗ₁ = Examples.bound_edges(s, ∂ₗ₀)
∂ᵣ₁ = Examples.bound_edges(s, ∂ᵣ₀)
∂ₗ₁₊ = Examples.adj_edges(s, ∂ₗ₀)
∂ᵣ₁₊ = Examples.adj_edges(s, ∂ᵣ₀)
∂ₑ₁₊ = Examples.adj_edges(s, ∂ₑ₀)
∂ₒ₁ = Examples.bound_edges(s, ∂ₒ₀)

∂ₗ₀₊ = unique(vcat(s[∂ₗ₁₊, :tgt], s[∂ₗ₁₊, :src]))
∂ᵣ₀₊ = unique(vcat(s[∂ᵣ₁₊, :tgt], s[∂ᵣ₁₊, :src]))
∂ₑ₀₊ = unique(vcat(s[∂ₑ₁₊, :tgt], s[∂ₑ₁₊, :src]))

# Heat diffusion constant in fluid
k₁ = 0.118
# Heat diffusion constant in cylinder
k₂ = 0.236980717779521 * 2

k_col = fill(k₁, ne(s))
k_col[cyl_edge] .= k₂
k = diagm(k_col)

density = 0.000210322
R₀ = 286.9
kᵥ = 0.085

# Pressure field for continuity
kₚ = 1600.0
kᵥₚ = 0.1

###################################
# BCs/Constants/Special Operators #
###################################

funcs = sym2func(sd)
funcs[:mcopy] = Dict(:operator => I(ne(sd)), :type => MatrixFunc())
funcs[:k] = Dict(:operator => k, :type => MatrixFunc())
funcs[:half] = Dict(:operator => 0.5 * I(nv(sd)), :type => MatrixFunc())
funcs[:third] = Dict(:operator => (1/3) * I(ne(sd)), :type => MatrixFunc())
funcs[:kₚ] = Dict(:operator => kₚ * I(nv(sd)), :type => MatrixFunc())
funcs[:R₀] = Dict(:operator => R₀ * I(nv(sd)), :type => MatrixFunc())
funcs[:kᵥ] = Dict(:operator => kᵥ * I(ne(sd)), :type => MatrixFunc())
funcs[:kᵥₚ] = Dict(:operator => kᵥₚ * I(nv(sd)), :type => MatrixFunc())
funcs[:dneg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:mult₀] = Dict(:operator => (x′,x,y) -> (x′ .= x .* y), :type => InPlaceFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₀] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())
funcs[:∂₀] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀₊] .= 288.15), :type => InPlaceFunc())
funcs[:∂₀₋] = Dict(:operator => (x′, x) -> (x′[∂ₗ₀] .= 0; x′[∂ᵣ₀] .= 0), :type => InPlaceFunc())
funcs[:∂ₗ₀₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₗ₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0; x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:∂₀ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₀] .= 0; x′[∂ₑ₀₊] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[cyl_edge] .= 0; x′[∂ₒ₁] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₀ₗ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀₊] .= 0), :type => InPlaceFunc())
funcs[:mask₀ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[cyl_inner] .= 0), :type => InPlaceFunc())

wedge_cache = init_wedge_ops(sd)
v2comp = comp_support(sd);
cache_mat = Dict(:t2c => tri2comp(s, v2comp), :e2c => edge2comp(s, v2comp), :cross => changes(sd, v2comp))
funcs[:∧₁₀′] = Dict(:operator => (x′, α, β) -> (x′ .= cp_2_1(β, α, cache_mat)), :type => InPlaceFunc())
funcs[:∧₁₁′] = Dict(:operator => (x′, α, β) -> (x′ .= pd_wedge(Val{(1,1)},s, α, β; wedge_cache...)), :type => InPlaceFunc())



##########################
# Set Initial Conditions #
##########################
c_objs = fill(288.15, nv(s))
c_objs[∂ₒ₀] .= 350.0

velocity(p) = [3.402,0.0,0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
v[cyl_edge] .= 0

p = [0.0 for p in s[:point]]
u0 = vcat(c_objs, v, p)

using JSON
# Importing initial smoothed velocity field
pre_cond = open("../final_cond_v2.json", "r") do f
    JSON.parse(f)
end
u0[(nv(s)+1):(nv(s)+ne(s))] .= pre_cond[(nv(s)+1):(nv(s)+ne(s))];


##################
# Run Simulation #
##################
# Compress current scheduling
new_dwd = deepcopy(exp_dwd)
Examples.contract_matrices!(new_dwd, funcs)

# Generate and run simulation
func, _ = gen_sim(new_dwd, funcs, sd; autodiff=false);
prob = ODEProblem(func, u0, (0.0, 5.0))
sol1 = solve(prob, Tsit5(), progress=true, progress_steps=100, saveat=0.05)


####################################
# Result Storage and Visualization #
####################################
using JSON
res = Dict{Float64, Vector}()
for i in 0.0:0.05:sol1.t[end]
    res[i] = sol1(i)
end
open("uncoupled_phys/sim_res2.json", "w") do f
    JSON.print(f, res)
end

using WriteVTK
for i in 101:(100+length(0.05:0.05:sol1.t[end]))
  tris = triangle_vertices(s)
  cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [tris[1][t],tris[2][t],tris[3][t] ]) for t in 1:ntriangles(s)];
  x =  [p[1] for p in s[:point]]
  y =  [p[2] for p in s[:point]]
  z =  [p[3] for p in s[:point]]
  vtkfile = vtk_grid("uncoupled_phys/su2_coupled_v1_$(lpad(i, 4, "0"))", x, y, z, cells)
  vel_form = ♯(sd, sol1(0.05*i - 5.0)[(1+nv(s)):(nv(s)+ne(s))], CombinatorialSpaces.DiscreteExteriorCalculus.PPSharp())

  vtkfile["temperature", VTKPointData()] = sol1(0.05*i - 5.0)[1:nv(s)]
  vtkfile["pressure", VTKPointData()] = sol1(0.05*i - 5.0)[(end-nv(s)+1):end]
  vtkfile["vel", VTKPointData()] = vel_form
  vtkfile["v_mag", VTKPointData()] = sqrt.(1e-7 .+ inv_hodge_star(Val{0}, sd)*pd_wedge(Val{(1,1)}, sd, sol1(0.05*i - 5.0)[(1+nv(s)):(nv(s)+ne(s))], ⋆(Val{1}, sd) * sol1(0.05*i - 5.0)[(1+nv(s)):(nv(s)+ne(s))]))
  vtk_save(vtkfile)
end

magnitudes(t) = sqrt.(1e-4 .+ inv_hodge_star(Val{0}, sd)*pd_wedge(Val{(1,1)}, sd, sol1(t)[(1+nv(s)):(nv(s)+ne(s))], ⋆(Val{1}, sd) * sol1(t)[(1+nv(s)):(nv(s)+ne(s))]))
times = range(0, sol1.t[end], length=150)
colors = [magnitudes(t) for t in times]
fig = Figure()
axis, scatter_thing = mesh(fig[1,1], s, color=colors[1], colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))

xlims!(axis, (-5,15))
ylims!(axis, (-5,5))
axis.aspect = AxisAspect(2.0)
Colorbar(fig[1,2], ob)
framerate = 30

record(fig, "vel_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
  scatter_thing.color = colors[i]
end
