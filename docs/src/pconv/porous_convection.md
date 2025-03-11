# Porous Convection

This Porous Convection model is based on the equations, constants and mesh structure given by ETH Zurich's course on Solving Partial Differential Equations in Parallel on GPUs [Lecture 4](https://pde-on-gpu.vaw.ethz.ch/lecture4).

```@setup INFO
include(joinpath(Base.@__DIR__, ".." , "..", "docinfo.jl"))
info = DocInfo.Info()
```

## Dependencies

```@example DEC
using CairoMakie
using Catlab
using CombinatorialSpaces
using ComponentArrays
using Decapodes
using DiagrammaticEquations
using Distributions
using GeometryBasics: Point2, Point3
using LinearAlgebra
using MLStyle
using OrdinaryDiffEq
using SparseArrays
using StaticArrays
using StatsBase

import CombinatorialSpaces.DiscreteExteriorCalculus: eval_constant_primal_form
```

## Porous Convection Decapode

The overall physics the model wishes to capture here are those of a fluid at different temperatures. Differences in the temperature of the fluid led to differences in the densities of that fluid which leads to convection of the fluid.

To drive the convection, the model has a cooling element at the top and a heating element at the bottom of our mesh. To avoid the need to determine the density of the fluid, the [Boussinesq approximation](https://en.wikipedia.org/wiki/Boussinesq_approximation_(buoyancy)) is used.

The model generates a divergence-free pressure field by solving  [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation) and the [Darcy flux](https://en.wikipedia.org/wiki/Darcy%27s_law) is a combination of the forces caused by pressure differences and buoyant forces. This Darcy flux leads to advection of the fluid and along with some heat diffusion, the overall change in temperature is captured.

Finer details of the physics can be found at the source listed above.

```@example DEC
Porous_Convection = @decapode begin
  (λ_ρ₀Cp, αρ₀, k_ηf, ϕ)::Constant
  (P, T, Adv, bound_T, bound_Ṫ)::Form0
  (g, qD)::Form1

  bound_T == adiabatic(T)

  # Darcy flux
  ρ == g ∧ (αρ₀ * bound_T)
  P == Δ⁻¹(δ(ρ))
  qD == -k_ηf * (d(P) - ρ)

  Adv == ⋆(interpolate(∧ᵈᵖ₁₁(⋆(d(bound_T)), qD)))
  Ṫ == -1/ϕ * Adv + λ_ρ₀Cp * Δ(bound_T)

  bound_Ṫ == tb_bc(Ṫ)

  ∂ₜ(T) == bound_Ṫ
end
infer_types!(Porous_Convection)
resolve_overloads!(Porous_Convection)
to_graphviz(Porous_Convection)
```

## The Mesh

Our mesh will be a triangulated grid mesh of width 40 units and height 20. We ues the circumcenter subdivision method as the triangulated grid's triangles are well-behaved and thus we can enjoy faster solve times.

```@example DEC
lx, ly = 40.0, 20.0
dx = dy = 0.4
s = triangulated_grid(lx, ly, dx, dy, Point3{Float64});
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point2{Float64}}(s);
subdivide_duals!(sd, Circumcenter());

fig = Figure() # hide
ax = CairoMakie.Axis(fig[1,1], aspect=2) # hide
wf = wireframe!(ax, s; linewidth=1) # hide
resize_to_layout!(fig) # hide
fig # hide
```

## Operators and Boundary Conditions

We set up our custom operators first. Since our system is fairly small, we can directly factorize our Laplacian to make solve the Poisson Equation fast and accurate.

As you may have noticed, all of our data is located on the primal mesh, including our temperature. While this means we can keep the Darcy flux on the primal elements, the Lie derivative needs to compute the advection is slightly complicated now. However, we can employ interpolation to faithfully map data from primal to dual triangle elements to sidestep this problem. For fine enough meshes, the error in doing so in negligible.

Along with our operators, we implement [adiabatic conditions](https://en.wikipedia.org/wiki/Adiabatic_process) on the side of our mesh as well as no-change conditions on the top/bottom of our mesh, where our cooling/heating elements reside. The adiabatic condition works by giving the boundary the values of the appropriate horizontal neighbor, which easy to do on the triangulated grid as its structure allows for easy lookup of this neighbor.

```@example DEC
Δ0 = Δ(0,sd);
fΔ0 = LinearAlgebra.factorize(Δ0);

mat = p2_d2_interpolation(sd)

left_wall_idxs = findall(p -> p[1] < dx, s[:point]);
right_wall_idxs = findall(p -> p[1] > lx - dx, s[:point]);

# For adiabatic conditions:
next_left_wall_idxs = left_wall_idxs .+ 1;
next_right_wall_idxs = right_wall_idxs .- 1;

# For no-change conditions
bottom_wall_idxs= findall(p -> p[2] == 0, s[:point]);
top_wall_idxs = findall(p -> p[2] == ly, s[:point]);

apply_tb_bc(x) = begin x[bottom_wall_idxs] .= 0; x[top_wall_idxs] .= 0; return x; end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ⁻¹ => x -> begin
      y = fΔ0 \ x
      # Constant changes in solution are valid
      y .-= minimum(y)
    end
    :adiabatic => x -> begin
      x[left_wall_idxs] .= x[next_left_wall_idxs]
      x[right_wall_idxs] .= x[next_right_wall_idxs]
      return x
    end
    :tb_bc => apply_tb_bc
    :interpolate => x -> mat * x
    _ => error("No operator $my_symbol found.")
  end
  return op
end

sim = eval(gensim(Porous_Convection))
f = sim(sd, generate, DiagonalHodge())
nothing # hide
```

## Initial Conditions

We set up a Gaussian heat disturbance at the center of our mesh to produce interesting convective behavior. We also set the cooling/heating elements.

```@example DEC
ΔT = 200.0

T_dist = MvNormal([lx/2.0, ly/2.0], [1/sqrt(2), 1/sqrt(2)])
T = [2 * ΔT * pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]
T[top_wall_idxs] .= -ΔT/2
T[bottom_wall_idxs] .= ΔT/2

fig = Figure() # hide
ax = CairoMakie.Axis(fig[1,1], aspect=2) # hide
msh = mesh!(ax, s; color = T, colormap=:jet, colorrange=(-ΔT/2, ΔT/2)) # hide
Colorbar(fig[1,2], msh) # hide
colsize!(fig.layout, 1, Aspect(1, 2.0)) # hide
resize_to_layout!(fig) # hide
fig # hide
```

Since we depend on gravity to drive the convection, as this is why lower densities rise against higher densities, we define gravity as a force moving in the downward y-direction and use `eval_constant_primal_form` to generate the appropriate 1-Form to represent this force.

Afterwards we establish a variety of physical constants.

```@example DEC
# Gravity
accl_g = 9.81
grav = SVector{3}([0.0, -accl_g, 0.0])
g = eval_constant_primal_form(sd, grav)

# Physical constants
Ra = 750
k_ηf = 1.0
αρ₀ = (1.0/accl_g)
ϕ = 0.1
λ_ρ₀Cp = 1/Ra*(accl_g*αρ₀*k_ηf*ΔT*ly/ϕ)

u₀ = ComponentArray(T=T, g=g)
constants = (k_ηf = k_ηf, αρ₀ = αρ₀, ϕ = ϕ, λ_ρ₀Cp = λ_ρ₀Cp)
```

## Solving the Porous Convection Equations

We simply plug in our initial conditions and generate simulation code and solve it using `Tsit5()`. Below is an output of the full simulation.

```@example DEC
tₑ = 0.7
prob = ODEProblem(f, u₀, (0, tₑ), constants)
soln = solve(prob, Tsit5(); saveat = 0.005)
soln.retcode
```

```@setup DEC
wdg10 = dec_wedge_product(Tuple{1, 0}, sd)
codif_1 = δ(1, sd)
d0 = dec_differential(0, sd)

hdg_1 = dec_hodge_star(1, sd)
dp_wdg_11 = dec_wedge_product_dp(Tuple{1,1}, sd)
inv_hdg_0 = dec_inv_hodge_star(0, sd)

function calculate_pressure(T, constants)
  fΔ0 \ (codif_1 * wdg10(g, constants.αρ₀ * T))
end

function calculate_advection(T, P, constants)
  darcy_flux = -constants.k_ηf * (d0 * P - wdg10(g, constants.αρ₀ * T))
  apply_tb_bc(-1/constants.ϕ *  inv_hdg_0 * mat * dp_wdg_11(hdg_1 * d0 * T, darcy_flux))
end

function calculate_diffusion(T, constants)
  apply_tb_bc(constants.λ_ρ₀Cp * Δ0 * T)
end

function compute_colorranges(length)
  values = hcat(map(range(0, soln.t[end], length=length)) do t
    T = soln(t).T
    P = calculate_pressure(T, constants)
    Adv = calculate_advection(T, P, constants)
    Diff = calculate_diffusion(T, constants)

    [minimum(P), maximum(P), minimum(Adv), maximum(Adv), minimum(Diff), maximum(Diff)]
  end...)

  percentile(values[1, :], 90), percentile(values[2, :], 90), # Pressure
  percentile(values[3, :], 90), percentile(values[4, :], 90), # Advection
  percentile(values[5, :], 75), percentile(values[6, :], 75) # Diffusion
end

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)
  P = @lift(calculate_pressure($T, constants))
  Adv = @lift(calculate_advection($T, $P, constants))
  Diff = @lift(calculate_diffusion($T, constants))

  colorranges = compute_colorranges(video_length)
  P_range = colorranges[1], colorranges[2]
  Adv_range = colorranges[3], colorranges[4]
  Diff_range = colorranges[5], colorranges[6]

  f = Figure()

  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Temperature at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s; color=T, colormap=:jet, colorrange=(-ΔT/2, ΔT/2))
  Colorbar(f[1,2], msh_T)

  ax_P = CairoMakie.Axis(f[2,1], title = @lift("Pressure at Time $(round($time, digits=3))"))
  msh_P = mesh!(ax_P, s; color=P, colormap=:jet, colorrange=P_range)
  Colorbar(f[2,2], msh_P)

  ax_Adv = CairoMakie.Axis(f[1,3], title = @lift("Advection at Time $(round($time, digits=3))"))
  msh_Adv = mesh!(ax_Adv, s; color=Adv, colormap=:jet, colorrange=Adv_range)
  Colorbar(f[1,4], msh_Adv)

  ax_Diff = CairoMakie.Axis(f[2,3], title = @lift("Diffusion at Time $(round($time, digits=3))"))
  msh_Diff = mesh!(ax_Diff, s; color=Diff, colormap=:jet, colorrange=Diff_range)
  Colorbar(f[2,4], msh_Diff)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end
```

```@example DEC
save_dynamics("Porous_Convection.mp4", 120)
```

!["Porous Convection Result"](Porous_Convection.mp4)

```@example INFO
DocInfo.get_report(info) # hide
```
