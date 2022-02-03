# Our developed libraries
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus


using Catlab.Present
using Catlab.Graphics
using Catlab.CategoricalAlgebra
using LinearAlgebra

using Decapodes.Simulations
using Decapodes.Examples
using Decapodes.Diagrams
using Decapodes.Schedules

# Julia community libraries
using MeshIO
using CairoMakie
using Decapodes.Debug
using DifferentialEquations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

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

	# Boundary Conditions (currently does not impact actual BCs)
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

# Extend to incompressible Navier Stokes

@present IncompressibleFlow2D <: DiffAdv2D begin
  P == V ⋅ ⋆₁(X)
  dV == (V ⊗ P) ⋅ L₁ ⋅ dneg₁ ⋅ ⋆₁⁻¹(X) ⋅ mask₁ₑ + V ⋅ Δ₁(X) ⋅ kᵥ + p ⋅ d₀(X) ⋅ neg₁
  V ⋅ ∂ₜ(Form1(X)) == dV
  dp == V ⋅ neg₁ ⋅ δ₁(X) ⋅ kₚ + (V ⊗ (p⋅⋆₀(X))) ⋅ L₀ ⋅ ⋆₀⁻¹(X) ⋅ mask₀ₗ ⋅ neg₀
  V ⋅ ∂₁ₛ == bc₁
  dC ⋅ ∂₀ₛ == bc₃
  dp ⋅ ∂₀₋ == bc₀
  dV ⋅ ∂₁ₗ₊ == bc₁
  p ⋅ ∂ₜ(Form0(X)) == dp
end

diag = eq_to_diagrams(IncompressibleFlow2D)
to_graphviz(diag ; edge_len="1.3")

dwd =diag2dwd(diag)
to_graphviz(dwd, orientation=LeftToRight)

exp_dwd = Examples.expand_dwd(dwd, gen_dec_rules())
to_graphviz(exp_dwd, orientation=LeftToRight)

# Load mesh

s = EmbeddedDeltaSet2D("../meshes/naca0012_8.stl")
sd = dual(s);

# Define non-default operators (multiplication by constants, boundary conditions)
k = 0.1
kₚ = 10.0
kᵥ = 1.0
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
funcs[:half] = Dict(:operator => 0.5 * I(ne(sd)), :type => MatrixFunc())
funcs[:kₚ] = Dict(:operator => kₚ * I(nv(sd)), :type => MatrixFunc())
funcs[:kᵥ] = Dict(:operator => kᵥ * I(ne(sd)), :type => MatrixFunc())
funcs[:dneg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())
funcs[:neg₀] = Dict(:operator => -1 * I(nv(sd)), :type => MatrixFunc())
funcs[:∂₀] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀] .= 0), :type => InPlaceFunc())
funcs[:∂₀₋] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀] .= 0; x′[∂ᵣ₀] .= 0), :type => InPlaceFunc())
funcs[:∂ₗ₀₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₗ₊] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0), :type => InPlaceFunc())
funcs[:∂₀ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₀] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₛ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₒ₁] .= 0), :type => InPlaceFunc())
funcs[:∂₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₁ₑ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₁₊] .= 0;  x′[∂ᵣ₁₊] .= 0), :type => InPlaceFunc())
funcs[:mask₀ₗ] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₗ₀₊] .= 0), :type => InPlaceFunc())

# Define initial conditions
c_objs = zeros(nv(s))
c_objs[∂ₒ₀] .= 1e-3
velocity(p) = [2.0 * exp(-(p[1]+50)/5),0.0,0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;
p = [0.0 for p in s[:point]]
u0 = vcat(c_objs, v, p)
GC.gc()

fig, ax, ob = wireframe(s)
ax.aspect = AxisAspect(3)
fig
save("mesh.svg", fig)

# Compress wiring diagram over contiguous matrix multiplications
cont_dwd = deepcopy(exp_dwd)
Examples.contract_matrices!(cont_dwd, funcs)

# Generate simulation function
func, _ = gen_sim(cont_dwd, funcs, sd; autodiff=false);

# Solve problem
prob = ODEProblem(func, u0, (0, 30))
sol = solve(prob, Tsit5(), progress=true, progress_steps=100, p=v);

# Example of debugging simulation

# Key for debugging simulation
sim_key(exp_dwd, orientation=LeftToRight)

exp_func, _ = gen_sim(exp_dwd, funcs, sd; autodiff=false);

# Show the flow of concentration as arrows
fig, ax, ob = draw_wire(s, sd, exp_dwd, exp_func, sol(30.0), 52; axisaspect=1, xlim=(40, 50), ylim=(-5, 5), use_arrows=true, n_arrows=1000)
fig
save("debug_conc_flow.svg")

# Generate animations of simulation

times = range(0, sol.t[end], length=150)
colors = [sol(t)#=[(end-nv(s)+1):end]=# for t in times]

figure, axis, scatter_thing = mesh(s, color=colors[1],
                                   colorrange=(0,1e-3))#minimum(vcat(colors...)),maximum(vcat(colors...))))
axis.aspect = AxisAspect(100.0/30.0)
framerate = 30

record(figure, "conc_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
  scatter_thing.color = colors[i]
end

times = range(0, sol.t[end], length=150)
colors = [sol(t)[(end-nv(s)+1):end] for t in times]

figure, axis, scatter_thing = mesh(s, color=colors[1],
                                   colorrange=(minimum(vcat(colors...)),maximum(vcat(colors...))))
axis.aspect = AxisAspect(100.0/30.0)
framerate = 30

record(figure, "press_flow.gif", collect(1:length(collect(times))); framerate = framerate) do i
  scatter_thing.color = colors[i]
end
