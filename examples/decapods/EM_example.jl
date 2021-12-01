# Our developed libraries
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus


using Catlab.Present
using Catlab.Graphics
using Catlab.CategoricalAlgebra
using LinearAlgebra

# Julia community libraries
using MeshIO
using CairoMakie
using DifferentialEquations
using Distributions
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Decapods.Simulations
using Decapods.Examples
using Decapods.Diagrams
using Decapods.Schedules
using Decapods.Debug

# Define used quantities
@present EM2DQuantities(FreeExtCalc2D) begin
    X::Space
    neg₁::Hom(Form1(X), Form1(X)) # -1
    E::Hom(munit(), Form1(X))     # electric field
    B::Hom(munit(), Form2(X))     # magnetic field
    c::Hom(Form2(X), Form2(X))   # μ₀ / ϵ₀ (scalar)
end

# Define Electromagnetic physics
@present EM2D <: EM2DQuantities begin
    B ⋅ ∂ₜ(Form2(X)) == E ⋅ neg₁ ⋅ d₁(X)
    E ⋅ ∂ₜ(Form1(X)) == B ⋅ c ⋅ ⋆₂(X) ⋅ dual_d₀(X) ⋅ ⋆₁⁻¹(X)
end

diag = eq_to_diagrams(EM2D)
to_graphviz(diag; edge_len = "1.3")
##

dwd = diag2dwd(diag)
to_graphviz(dwd, orientation = LeftToRight)
##

s = EmbeddedDeltaSet2D("../meshes/pipe_fine.stl")
sd = dual(s);

c = -4.0
# c = -1 * (1.2566 / 8.8542)
funcs = sym2func(sd)
funcs[:c] = Dict(:operator => c * I(ntriangles(s)), :type => MatrixFunc())
funcs[:neg₁] = Dict(:operator => -1 * I(ne(sd)), :type => MatrixFunc())

# Initial Conditions
bField(x) = begin
    mag = (x[1] < -30) ? sin(x[1] / 3.0) : 0.0
    mag
end

E = zeros(Float64, ne(s))
B = [bField(sd[sd[t, :tri_center], :dual_point]) for t in 1:ntriangles(s)]

# Generate leapfrog simulation
f1, f2, dwd1, dwd2 = gen_leapfrog(diag, funcs, sd, [:B, :E]);
tspan = (0.0,10.0)
dyn_prob = DynamicalODEProblem(f1,f2,B,E,tspan)
sol = solve(dyn_prob, VerletLeapfrog(); dt=0.01, progress = true, progress_steps=100);

exp_func, _ = gen_sim(dwd, funcs, sd; autodiff = false);

# 
fig, ax, ob = draw_wire(s, sd, dwd, exp_func, sol[end-1], 6)
fig

# Plot solution
B_range = 1:ntriangles(s)

fig, ax, ob = draw_wire(s, sd, dwd, exp_func, sol(10.0), 3)#; colorrange = (-1e-0, 1e-0))
ax.aspect=AxisAspect(3.0)
fig

# Record solution
times = range(1e-4, sol.t[end], length = 300)
colors = [vcat([[v, v, v] for v in sol(t)[B_range]]...) for t in times]
colorrange = maximum(vcat(colors...))

framerate = 30

record(fig, "magnetic_field.gif", collect(1:length(collect(times))); framerate = framerate) do i
    ob.color = colors[i]
    ob.colorrange= (-4, 4)
end
