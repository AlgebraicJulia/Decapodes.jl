using Catlab, Catlab.CategoricalAlgebra
using CombinatorialSpaces, CombinatorialSpaces.SimplicialSets
using GeometryBasics: Point3

# TODO: Is there an easier way to ensure compatability with Makie plotting than
# adding this dependency now just to use their Point type?
Point3D = Point3{Float64}

"""
    makeSphere(minLat, maxLat, dLat, minLong, maxLong, dLong, radius)

Construct a spherical mesh (inclusively) bounded by the given latitudes and
longitudes, discretized at dLat and dLong intervals, at the given radius from
Earth's center.

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
# Regular octahedron.
julia> s, npi, spi = makeSphere(0, 180, 90, 0, 360, 90, 1)
````
```julia-repl
# 72 points along the unit circle on the x-y plane.
julia> s, npi, spi = makeSphere(90, 90, 0, 0, 360, 5, 1)
````
```julia-repl
# 72 points along the equator at 0km from Earth's surface.
julia> s, npi, spi = makeSphere(90, 90, 1, 0, 360, 5, 6371)
````
```julia-repl
# TIE-GCM grid at 90km altitude (with no poles,   i.e. a bulbous cylinder).
julia> s, npi, spi = makeSphere(5, 175, 5, 0, 360, 5, 6371+90)
````
```julia-repl
# TIE-GCM grid at 90km altitude (with South pole, i.e. a bowl).
julia> s, npi, spi = makeSphere(5, 180, 5, 0, 360, 5, 6371+90)
````
```julia-repl
# TIE-GCM grid at 90km altitude (with poles,      i.e. a sphere).
julia> s, npi, spi = makeSphere(0, 180, 5, 0, 360, 5, 6371+90)
````
```julia-repl
# The Northern hemisphere of the TIE-GCM grid at 90km altitude.
julia> s, npi, spi = makeSphere(0, 180, 5, 0, 360, 5, 6371+90)
````
"""
function makeSphere(minLat, maxLat, dLat, minLong, maxLong, dLong, radius)
  if (   !(0 ≤ minLat ≤ 180)  || !(0 ≤ maxLat ≤ 180)
      || !(0 ≤ minLong ≤ 360) || !(0 ≤ minLong ≤ 360)
      ||  (maxLat < minLat)   ||  (maxLong < minLong))
    throw(ArgumentError("Mins must be less than Maxs, lats must be in [0,180],"*
                        " and longs must be in [0,360]."))
  end
  sph2car(ρ,ϕ,θ) = (ρ*sind(θ)*cosd(ϕ),
                    ρ*sind(θ)*sind(ϕ),
                    ρ*cosd(θ))
  if (minLat == maxLat && dLat == 0)
    dLat = 1 # User wants a just one parallel. a:0:a is not valid julia.
  end
  s = EmbeddedDeltaSet2D{Bool, Point3D}()
  # Neither pole counts as a Meridian.
  connect_north_pole = false
  connect_south_pole = false
  if (minLat == 0)
    # Don't create num_meridians points at the North pole.
    minLat += dLat
    connect_north_pole = true
  end
  if (maxLat == 180)
    # Don't create num_meridians points at the South pole.
    maxLat -= dLat
    connect_south_pole = true
  end
  connect_long = false
  if (maxLong == 360)
    maxLong -= dLong
    connect_long = true
  end
  # TODO: Should we warn the user if the stitching edges are shorter than the
  # rest?
  num_parallels = length(minLat:dLat:maxLat)
  num_meridians = length(minLong:dLong:maxLong)
  ρ = radius
  # Add points one parallel at a time.
  for θ in minLat:dLat:maxLat
    vertex_offset = nv(s)+1
    add_vertices!(s, num_meridians,
                  point=map(minLong:dLong:maxLong) do ϕ
                    Point3D(sph2car(ρ,ϕ,θ)...)
                  end)
    # Connect this parallel.
    if (connect_long)
      add_sorted_edge!(s, vertex_offset+num_meridians-1, vertex_offset)
    end
    # Don't make triangles with the previous parallel if there isn't one.
    if (vertex_offset == 1)
      add_sorted_edges!(s,
                        vertex_offset:vertex_offset+num_meridians-2,
                        vertex_offset+1:vertex_offset+num_meridians-1)
      continue
    end
    # Add the triangles.
    foreach(vertex_offset  :vertex_offset+num_meridians-2,
            vertex_offset+1:vertex_offset+num_meridians-1,
            vertex_offset-num_meridians:vertex_offset-2) do i,j,k
      glue_sorted_triangle!(s, i, j, k)
    end
    foreach(vertex_offset+1:vertex_offset+num_meridians-1,
            vertex_offset-num_meridians:vertex_offset-2,
            vertex_offset-num_meridians+1:vertex_offset-1) do i,j,k
      glue_sorted_triangle!(s, i, j, k)
    end
    # Connect with the previous parallel.
    if (connect_long)
      glue_sorted_triangle!(s, vertex_offset+num_meridians-1,
                            vertex_offset, vertex_offset-1)
      glue_sorted_triangle!(s, vertex_offset-num_meridians,
                            vertex_offset, vertex_offset-1)
    end
  end
  # Add the North and South poles.
  north_pole_idx = 0
  if (connect_north_pole)
    ϕ, θ = 0, 0
    add_vertex!(s, point=Point3D(sph2car(ρ,ϕ,θ)...))
    north_pole_idx = nv(s)
    foreach(1:num_meridians-1, 2:num_meridians) do i,j
      glue_sorted_triangle!(s, north_pole_idx, i, j)
    end
    if (connect_long)
      glue_sorted_triangle!(s, north_pole_idx, 1, num_meridians)
    end
  end
  south_pole_idx = 0
  if (connect_south_pole)
    south_parallel_start = num_meridians*(num_parallels-1)+1
    ϕ, θ = 0, 180
    add_vertex!(s, point=Point3D(sph2car(ρ,ϕ,θ)...))
    south_pole_idx = nv(s)
    foreach(south_parallel_start:south_parallel_start+num_meridians-2,
            south_parallel_start+1:south_parallel_start+num_meridians-1) do i,j
      glue_sorted_triangle!(s, south_pole_idx, i, j)
    end
    if (connect_long)
      glue_sorted_triangle!(s, south_pole_idx, south_parallel_start,
                            south_parallel_start+num_meridians-1)
    end
  end
  return s, north_pole_idx, south_pole_idx
end

## Tests
#using Test
#ρ = 6371+90
#s, npi, spi = makeSphere(0, 180, 5, 0, 360, 5, ρ)
#magnitude = (sqrt ∘ (x -> foldl(+, x*x)))
#
## The definition of a discretization of a sphere of unspecified radius.
#ρ′ = magnitude(s[:point][begin])
#@test all(isapprox.(magnitude.(s[:point]), ρ′))
#
## The definition of a discretization of a sphere of radius ρ.
#@test all(isapprox.(magnitude.(s[:point]), ρ))
#
## Some properties of a regular octahedron.
#◀▶, npi, spi = makeSphere(0, 180, 90, 0, 360, 90, 1)
#@test nv(◀▶) == 6
#@test ne(◀▶) == 12
#@test length(triangles(◀▶)) == 8

