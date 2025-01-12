using Catlab
using CombinatorialSpaces
using DiagrammaticEquations
using Decapodes
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
import CairoMakie: wireframe, mesh, Figure, Axis
using ComponentArrays

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

# We use the model equations as stated here:
# https://github.com/JuliaParallel/julia-hpc-tutorial-sc24/blob/main/parts/gpu/gray-scott.ipynb
# Initial conditions were based off those given here:
# https://itp.uni-frankfurt.de/~gros/StudentProjects/Projects_2020/projekt_schulz_kaefer/#header
GrayScott = @decapode begin
  (U, V)::Form0
  (UV2)::Form0
  (U̇, V̇)::Form0
  (f, k, rᵤ, rᵥ)::Constant
  B::Constant

  UV2 == (U .* (V .* V))
  lap_U == mask(Δ(U), B)
  lap_V == mask(Δ(V), B)

  U̇ == rᵤ * lap_U - UV2 + f * (1 .- U)
  V̇ == rᵥ * lap_V + UV2 - (f + k) .* V
  ∂ₜ(U) == U̇
  ∂ₜ(V) == V̇
end

n = 100
h = 1

s = triangulated_grid(n,n,h,h,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s);
subdivide_duals!(sd, Circumcenter());

sim = eval(gensim(GrayScott))

left_wall_idxs = findall(x -> x[1] <= h, s[:point])
right_wall_idxs = findall(x -> x[1] >= n - h, s[:point])
top_wall_idxs = findall(y -> y[2] == 0.0, s[:point])
bot_wall_idxs = findall(y -> y[2] == n, s[:point])

wall_idxs = unique(vcat(left_wall_idxs, right_wall_idxs, top_wall_idxs, bot_wall_idxs))
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :mask => (x,y) -> begin
      x[wall_idxs] .= y
      x
    end
    _ => error("Unmatched operator $my_symbol")
  end
end

fₘ = sim(sd, generate, DiagonalHodge())

init_multi = 0.5

U = rand(0.0:0.001:0.1, nv(sd))
V = zeros(nv(sd))

mid = div(n, 2)

mid_p = Point2D(mid, mid)

init = map(p -> if norm(p - mid_p, Inf) <= 5; 1.0 .* init_multi; else 0.0; end, sd[:point])

# Set up an initial small disturbance
U .+= init
V .+= 0.5 * init

u₀ = ComponentArray(U=U,V=V)

f = 0.055
k = 0.062
constants_and_parameters = (
  rᵤ = 0.16,
  rᵥ = 0.08,
  f = f,
  k = k,
  B = 0)

# fig = Figure();
# ax = CairoMakie.Axis(fig[1,1], aspect=1, title = "Initial value of U")
# msh = CairoMakie.mesh!(ax, s, color=U, colormap=:jet, colorrange=(extrema(U)))
# Colorbar(fig[1,2], msh)
# display(fig)

# fig = Figure()
# ax = CairoMakie.Axis(fig[1,1], aspect=1, title = "Initial value of V") # hide
# msh = CairoMakie.mesh!(ax, s, color=V, colormap=:jet, colorrange=extrema(V)) # hide
# Colorbar(fig[1,2], msh)
# fig

tₑ = 10_000

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

function save_dynamics(save_file_name)
  time = Observable(0.0)
  u = @lift(soln($time).U)
  f = Figure()
  ax_U = CairoMakie.Axis(f[1,1], title = @lift("Concentration of U at Time $($time)"))

  msh_U = mesh!(ax_U, s, color=u, colormap=:jet, colorrange=(0, 1.1))
  Colorbar(f[1,2], msh_U)

  timestamps = range(0, tₑ, step=50)
  record(f, save_file_name, timestamps; framerate = 30) do t
    time[] = t
  end
end

save_dynamics("gs_f=$(f)_k=$(k).mp4")