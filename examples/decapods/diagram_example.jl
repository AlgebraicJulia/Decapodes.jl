# Our developed libraries
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus

using Decapods.Simulations
using Decapods.Examples
using Decapods.Diagrams
using Decapods.Schedules

using Catlab.Present
using Catlab.Graphics

# Julia community libraries
using MeshIO
using CairoMakie
using DifferentialEquations
using LinearAlgebra

@present Diffusion2DQuantities(FreeExtCalc2D) begin
  X::Space
  C::Hom(munit(), Form0(X))     # concentration
  ϕ::Hom(munit(), DualForm1(X)) # negative diffusion flux
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  ∂₀::Hom(Form0(X), Form0(X))
  zero₀::Hom(munit(), Form0(X))
end

@present Diffusion2D <: Diffusion2DQuantities begin
  # Fick's first law
  ϕ == C ⋅ d₀(X) ⋅ k ⋅ ⋆₁(X)
  # Diffusion equation
  C ⋅ ∂ₜ(Form0(X)) == ϕ ⋅ dual_d₁(X) ⋅ ⋆₀⁻¹(X)
  C ⋅ ∂₀ == zero₀
end;

@present DiffAdv2DQuantities(FreeExtCalc2D) begin
  X::Space
  C::Hom(munit(), Form0(X))     # concentration
  V::Hom(munit(), Form1(X))     # Flow field
  ϕ::Hom(munit(), DualForm1(X)) # negative diffusion flux
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  ∂₀::Hom(Form0(X), Form0(X))
  zero₀::Hom(munit(), Form0(X))
end

@present DiffAdv2D <: DiffAdv2DQuantities begin
  # Fick's first law
  ϕ == C ⋅ d₀(X) ⋅ k ⋅ ⋆₁(X)
  # Diffusion equation
  C ⋅ ∂ₜ(Form0(X)) == ϕ ⋅ dual_d₁(X) ⋅ ⋆₀⁻¹(X) + (V⊗(C⋅⋆₀(X))) ⋅ L₀ ⋅ ⋆₀⁻¹(X)
  C ⋅ ∂₀ == zero₀
end

diag = eq_to_diagrams(Diffusion2D)
to_graphviz(diag ; edge_len="1.3")

new_dwd = diag2dwd(diag)
to_graphviz(new_dwd, orientation=LeftToRight)

s = EmbeddedDeltaSet2D("../meshes/naca0012_8.stl")
sd = dual(s);
fig, ax, ob = wireframe(s)
ax.aspect = AxisAspect(3)
fig

funcs = sym2func(sd)
const k = 0.3

∂₀ = Examples.boundary_inds(Val{0}, s)
∂ₗ₀ = ∂₀[findall(p-> p[1] <= -50.0, s[∂₀, :point])]
∂ᵣ₀ = ∂₀[findall(p-> p[1] >= 50.0, s[∂₀, :point])]
∂ₜ₀ = ∂₀[findall(p-> p[2] >= 15.0, s[∂₀, :point])]
∂ᵦ₀ = ∂₀[findall(p-> p[2] <= -15.0, s[∂₀, :point])]
∂ₑ₀ = vcat(∂ₗ₀, ∂ᵣ₀, ∂ₜ₀, ∂ᵦ₀)

funcs[:k] = Dict(:operator => k * I(ne(sd)), :type => MatrixFunc())
funcs[:∂₀] = Dict(:operator => (x′,x) -> (x′ .= x; x′[∂ₑ₀] .= 0), :type => InPlaceFunc())

diffusion_func, _ = gen_sim(new_dwd, funcs, sd; params=[:V]);

using Distributions
#Gaussian curve
x1 = -30.0
x2 = -20.0
c_dist = MvNormal([(x1+x2)/2, 0.0], [4.0, 4.0])
c_gauss = [pdf(c_dist, [p[1], p[2]]) for p in s[:point]]
max_gauss = maximum(c_gauss)
c_gauss ./= max_gauss


u0 = c_gauss
velocity(p) = [-3.0, 0.0, 0.0] #[p[2],-p[1],0.0]
v = ♭(sd, DualVectorField(velocity.(sd[triangle_center(sd),:dual_point]))).data;

prob = ODEProblem(diffusion_func, u0, (0, 10.0), v)
sol = solve(prob, Tsit5());

fig, ax, ob = mesh(s, color = sol(10), colorrange=(0,1))
ax.aspect = AxisAspect(3)
fig

times = range(0,10.0, length=150)
colors = [sol(t) for t in times]

figure, axis, scatter_thing = mesh(s, color=colors[1],
                                   colorrange=(minimum(colors[1]),maximum(colors[1])))
axis.aspect = AxisAspect(100.0/30.0)
framerate = 30

record(figure, "flow_conc.gif", collect(1:length(collect(times))); framerate = framerate) do i
  scatter_thing.color = colors[i]
end
