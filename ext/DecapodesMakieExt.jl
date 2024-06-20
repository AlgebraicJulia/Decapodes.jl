module DecapodesMakieExt

using CombinatorialSpaces
using Makie
import Decapodes: record_dynamics

"""
    record_dynamics(dynamics, mesh::HasDeltaSet, path::AbstractString, iter; use_colorbar=true, figure_args=NamedTuple(), axis_args=NamedTuple(), mesh_args=NamedTuple(), colorbar_args=NamedTuple(), record_args=NamedTuple())

The Decapodes implementation of Makie's `record` function. This function is meant serve as a helper function
to quickly generate animations based on the resulting dynamics of a Decapodes.

**Arguments:**

`dynamics`: A function that accepts a time and returns the value of a particular state variable at that time.

`mesh`: A CombinatorialSpaces mesh upon which the Decapode was run.

`path`: The path of the created file recording.

`iter`: The range of times where `dynamics` will be sampled.

**Keyword arguments:**

`use_colorbar = true`: Determines whether display a colorbar alongside the recording.

`figure_args`: Keyword arguments passed into Makie's `Figure` funcion.

`axis_args`: Keyword arguments passed into Makie's `Axis` funcion.

`mesh_args`: Keyword argumentsp assed into Makie's `mesh!` function. Note that `color` is already used to take the value of dynamics at time `t`.

`colorbar_args`: Keyword arguments passed into Makie's `Colorbar` funcion. Note that if `use_colorbar = false` these will have no effect.

`record_args`: Keyword arguments passed into Makie's `record` funcion.

"""
function record_dynamics(dynamics, mesh::HasDeltaSet, path::AbstractString, iter; use_colorbar=true, figure_args=NamedTuple(), axis_args=NamedTuple(), mesh_args=NamedTuple(), colorbar_args=NamedTuple(), record_args=NamedTuple())

  time = Observable(0.0)
  lift_dynam = @lift(dynamics($time))

  figure = Makie.Figure(;figure_args...)
  ax = Makie.Axis(figure[1,1]; axis_args...)

  makie_msh = Makie.mesh!(ax, mesh; color=lift_dynam, mesh_args...)

  use_colorbar ? Makie.Colorbar(figure[1,2], makie_msh; colorbar_args...) : nothing

  Makie.record(figure, path, iter; record_args...) do t
    time[] = t
  end
end

end
