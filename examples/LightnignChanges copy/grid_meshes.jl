using CombinatorialSpaces
using GeometryBasics: Point2, Point3
using LinearAlgebra: diagm
using GLMakie
Point2D = Point2{Float64}
Point3D = Point3{Float64}

# Note that dx will be slightly less than what is given, since points are
# compressed to lie within max_x.
function triangulated_grid(max_x, max_y, dx, dy)

  s = EmbeddedDeltaSet2D{Bool, Point3D}()

  # Place equally-spaced points in a max_x by max_y rectangle.
  coords = map(x -> Point3D(x..., 0), Iterators.product(0:dx:max_x, 0:dy:max_y))
  # Perturb every other row right by half a dx.
  map!(coords[:, 2:2:end], coords[:, 2:2:end]) do row
    row .+ [dx/2, 0, 0]
  end
  # The perturbation moved the right-most points past max_x, so compress along x.
  map!(coords, coords) do coord
    diagm([max_x/(max_x+dx/2), 1, 1]) * coord
  end

  add_vertices!(s, length(coords), point = vec(coords))

  nx = length(0:dx:max_x)

  # Matrix that stores indices of points.
  idcs = reshape(eachindex(coords), size(coords))
  # Only grab vertices that will be the bottom-left corner of a subdivided square.
  idcs = idcs[begin:end-1, begin:end-1]

  # Subdivide every other row along the opposite diagonal.
  for i in idcs[:, begin+1:2:end]
    glue_sorted_triangle!(s, i, i+nx, i+nx+1)
    glue_sorted_triangle!(s, i, i+1, i+nx+1)
  end
  for i in idcs[:, begin:2:end]
    glue_sorted_triangle!(s, i, i+1, i+nx)
    glue_sorted_triangle!(s, i+1, i+nx, i+nx+1)
  end

  # Orient and return.
  s[:edge_orientation]=true
  orient!(s)
  s
end
s = triangulated_grid(10, 10, 1, 1)
GLMakie.wireframe(s)
s = triangulated_grid(400, 100, 10, 5)
GLMakie.wireframe(s)

volumes = [CombinatorialSpaces.volume(2, s, i) for i in triangles(s)]
allequal(volumes)
sum(volumes) ≈ max_x*max_y

# Create a triangulated grid.
for i in eachindex(0:dy:max_y)[begin:end-1]
  for j in eachindex(0:dx:max_x)[begin:end-1]
    if i % 2 == 1
      glue_sorted_triangle!(s, (i-1)*(nx)+j, i*(nx)+j, i*(nx)+j+1)
      glue_sorted_triangle!(s, (i-1)*(nx)+j, (i-1)*(nx)+j+1, i*(nx)+j+1)
    else
      glue_sorted_triangle!(s, (i-1)*(nx)+j, (i-1)*(nx)+j+1, i*(nx)+j)
      glue_sorted_triangle!(s, (i-1)*(nx)+j+1, i*(nx)+j, i*(nx)+j+1)
    end
  end
end

s[:edge_orientation]=true
orient!(s)
s

GLMakie.wireframe(s)
#GLMakie.mesh(s, color=ones(nv(s)))
