using Catlab, Catlab.CategoricalAlgebra
using CombinatorialSpaces, CombinatorialSpaces.SimplicialSets
using Makie: Point3

Point3D = Point3{Float64}

"""
    makeSphere(minLat, maxLat, dLat, minLong, maxLong, dLong, radius,
      connectLong)

Construct a spherical mesh (inclusively) bounded by the given latitudes and
longitudes, discretized at dLat and dLong intervals, at the given radius from
Earth's center. If connectLong is true, then there are edges from the points at
minLong to the points at maxLong.

We say that:
- 90°N is 0
- 90°S is 180
- Prime Meridian is 0
- 10°W is 355

We say that:
- (x=0,y=0,z=0) is at the center of the sphere
- the x-axis points toward 0°,0°
- the y-axis points toward 90°E,0° TODO: Double-check this.
- the z-axis points toward the North Pole

# References:
[List of common coordinate transformations](https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations?oldformat=true#From_spherical_coordinates)

# Examples
```julia-repl
# 180 points along the unit circle on the x-y plane.
julia> s = makeSphere(90, 90, 1, 0, 175, 5, 1, true)
````
```julia-repl
# 180 points along the equator at 0km from Earth's surface.
julia> s = makeSphere(90, 90, 1, 0, 175, 5, 6371, true)
````
```julia-repl
# TIE-GCM grid
julia> s = makeSphere(5, 175, 5, 0, 175, 5, 6371+90, true)
````
"""
function makeSphere(minLat, maxLat, dLat, minLong, maxLong, dLong, radius,
    connectLong)
  s = EmbeddedDeltaSet2D{Bool, Point3D}()
  # Add points one parallel at a time.
  for θ in minLat:dLat:maxLat
    ρ = radius
    add_vertices!(s, length(minLong:dLong:maxLong),
                  point=map(minLong:dLong:maxLong) do ϕ
                    Point3D(ρ*sind(θ)*cosd(ϕ),
                            ρ*sind(θ)*sind(ϕ),
                            ρ*cosd(θ))
                  end)
  end
  return s
end

