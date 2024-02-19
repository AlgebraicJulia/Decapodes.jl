# This is a helper script that generates a triangulated grid.

using CombinatorialSpaces
using GeometryBasics: Point2, Point3
using LinearAlgebra: diagm
Point2D = Point2{Float64}
Point3D = Point3{Float64}

# Note that dx will be slightly less than what is given, since points are
# compressed to lie within max_x.
function triangulated_grid(max_x, max_y, dx, dy, point_type)

  s = EmbeddedDeltaSet2D{Bool, point_type}()

  ## Place equally-spaced points in a max_x by max_y rectangle.
  coords = point_type == Point3D ? map(x -> point_type(x..., 0), Iterators.product(0:dx:max_x, 0:dy:max_y)) : map(x -> point_type(x...), Iterators.product(0:dx:max_x, 0:dy:max_y))
  ## Perturb every other row right by half a dx.
  coords[:, 2:2:end] = map(coords[:, 2:2:end]) do row
    if point_type == Point3D
      row .+ [dx/2, 0, 0]
    else
      row .+ [dx/2, 0]
    end
  end
  ## The perturbation moved the right-most points past max_x, so compress along x.
  map!(coords, coords) do coord
    if point_type == Point3D
      diagm([max_x/(max_x+dx/2), 1, 1]) * coord
    else
      diagm([max_x/(max_x+dx/2), 1]) * coord
    end
  end

  add_vertices!(s, length(coords), point = vec(coords))

  nx = length(0:dx:max_x)

  ## Matrix that stores indices of points.
  idcs = reshape(eachindex(coords), size(coords))
  ## Only grab vertices that will be the bottom-left corner of a subdivided square.
  idcs = idcs[begin:end-1, begin:end-1]
  
  ## Subdivide every other row along the opposite diagonal.
  for i in idcs[:, begin+1:2:end]
    glue_sorted_triangle!(s, i, i+nx, i+nx+1)
    glue_sorted_triangle!(s, i, i+1, i+nx+1)
  end
  for i in idcs[:, begin:2:end]
    glue_sorted_triangle!(s, i, i+1, i+nx)
    glue_sorted_triangle!(s, i+1, i+nx, i+nx+1)
  end

  ## Orient and return.
  s[:edge_orientation]=true
  orient!(s)
  s
end
