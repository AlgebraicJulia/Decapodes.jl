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
- the y-axis points toward 90°E,0°
- the z-axis points toward the North Pole

# References:
[List of common coordinate transformations](https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations?oldformat=true#From_spherical_coordinates)

# Examples
```julia-repl
# 72 points along the unit circle on the x-y plane.
julia> s = makeSphere(90, 90, 0, 0, 355, 5, 1, true)
````
```julia-repl
# 72 points along the equator at 0km from Earth's surface.
julia> s = makeSphere(90, 90, 1, 0, 355, 5, 6371, true)
````
```julia-repl
# TIE-GCM grid at 90km altitude.
julia> s = makeSphere(5, 175, 5, 0, 355, 5, 6371+90, true)
````
"""
function makeSphere(minLat, maxLat, dLat, minLong, maxLong, dLong, radius,
    connectLong)
  if (   !(0 ≤ minLat ≤ 180)  || !(0 ≤ maxLat ≤ 180)
      || !(0 ≤ minLong ≤ 360) || !(0 ≤ minLong ≤ 360)
      ||  (maxLat < minLat)   ||  (maxLong < minLong))
    throw(ArgumentError(""))
  end
  if (minLat == maxLat && dLat == 0)
    dLat = 1 # User wants a unit circle. a:0:a is not valid julia.
  end
  s = EmbeddedDeltaSet2D{Bool, Point3D}()
  # Add points one parallel at a time.
  numMeridians = length(minLong:dLong:maxLong)
  for θ in minLat:dLat:maxLat
    vertexOffset = nv(s)+1
    ρ = radius
    add_vertices!(s, numMeridians,
                  point=map(minLong:dLong:maxLong) do ϕ
                    Point3D(ρ*sind(θ)*cosd(ϕ),
                            ρ*sind(θ)*sind(ϕ),
                            ρ*cosd(θ))
                  end)
    # Connect this parallel.
    if (connectLong)
      add_sorted_edge!(s, vertexOffset+numMeridians-1, vertexOffset)
    end
    # Don't make triangles with the previous parallel if there isn't one.
    if (vertexOffset == 1)
      add_sorted_edges!(s,
                        vertexOffset:vertexOffset+numMeridians-2,
                        vertexOffset+1:vertexOffset+numMeridians-1)
      continue
    end
    # Add the triangles.
    foreach(vertexOffset  :vertexOffset+numMeridians-2,
            vertexOffset+1:vertexOffset+numMeridians-1,
            vertexOffset-numMeridians:vertexOffset-2) do i,j,k
      glue_sorted_triangle!(s, i, j, k)
    end
    foreach(vertexOffset+1:vertexOffset+numMeridians-1,
            vertexOffset-numMeridians:vertexOffset-2,
            vertexOffset-numMeridians+1:vertexOffset-1) do i,j,k
      glue_sorted_triangle!(s, i, j, k)
    end
    # Connect with the previous parallel.
    if (connectLong)
      glue_sorted_triangle!(s, vertexOffset+numMeridians-1,
                            vertexOffset, vertexOffset-1)
      glue_sorted_triangle!(s, vertexOffset-numMeridians,
                            vertexOffset, vertexOffset-1)
    end
  end
  return s
end

