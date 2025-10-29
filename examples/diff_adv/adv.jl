using Decapodes
using DiagrammaticEquations
using CombinatorialSpaces
using GeometryBasics
using MLStyle
using ComponentArrays
using OrdinaryDiffEq
using CairoMakie
using Distributions

INITCOND = "test"
RESOLUTION = 0.5

lx = 100
ly = 30

s = triangulated_grid(lx, ly, RESOLUTION, RESOLUTION, Point3D)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

Adv = @decapode begin
    T::Form0
    u::Form1
    ∂ₜ(T) == -(⋆(d(⋆(T ∧ u))))
end
infer_types!(Adv)
resolve_overloads!(Adv)

simulate = evalsim(Adv)
fₘ = simulate(sd, nothing)

if INITCOND == "square"
    T = map(sd[:point]) do p
        if 40 <= p[1] <= 60 && 40 <= p[2] <= 60
            return 1
        else
            return 0
        end
    end
else
    T_dist = MvNormal([lx/2, ly/2], [5, 5])
    T = [pdf(T_dist, [p[1], p[2]]) for p in sd[:point]]
end
T₀ = T

fig = Figure();
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=T, colormap=:jet, colorrange=extrema(T))
Colorbar(fig[1,2], msh)
fig

u = eval_constant_primal_form(sd, Point3d(1, 0, 0))

u₀ = ComponentArray(T=T,u=u)

constants_and_parameters = ()

tₑ = 20

@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@info("Done")

FILENAME = "$(INITCOND)_mesh_$(RESOLUTION)_res"

function save_dynamics(save_file_name, video_length = 30)
  time = Observable(0.0)

  T = @lift(soln($time).T)
  f = Figure()
  ax_T = CairoMakie.Axis(f[1,1], title = @lift("Temperature at Time $(round($time, digits=3))"))
  msh_T = mesh!(ax_T, s; color=T, colormap=:jet, colorrange = extrema(T₀))
  Colorbar(f[1,2], msh_T)

  timestamps = range(0, soln.t[end], length=video_length)
  record(f, save_file_name, timestamps; framerate = 15) do t
    time[] = t
  end
end

save_dynamics("$(FILENAME).mp4")

# fig = Figure();
# ax = CairoMakie.Axis(fig[1,1])
# msh = CairoMakie.mesh!(ax, s, color=soln.u[end].T, colormap=:jet, colorrange=extrema(T))
# Colorbar(fig[1,2], msh)
# fig

energy_conv = map(t -> sum(soln(t).T), soln.t) .- sum(T)

fig = Figure(size=(2000, 1000));
ax = CairoMakie.Axis(fig[1,1])
msh = CairoMakie.mesh!(ax, s, color=soln.u[end].T, colormap=:jet, colorrange=extrema(T))
Colorbar(fig[1,2], msh)
ax2 = CairoMakie.Axis(fig[1,3]; title = "Energy conservation with $(INITCOND) mesh")
CairoMakie.plot!(ax2, soln.t, energy_conv)
save("$(FILENAME).png", fig)