# These are helper functions that return the boundary indices of a given mesh.

function boundary_inds(::Type{Val{1}}, s)
  collect(findall(x -> x != 0, boundary(Val{2},s) * fill(1,ntriangles(s))))
end

function boundary_inds(::Type{Val{0}}, s)
  ∂1_inds = boundary_inds(Val{1}, s)
  unique(vcat(s[∂1_inds,:∂v0],s[∂1_inds,:∂v1]))
end

function boundary_inds(::Type{Val{2}}, s)
  ∂1_inds = boundary_inds(Val{1}, s)
  inds = map([:∂e0, :∂e1, :∂e2]) do esym
    vcat(incident(s, ∂1_inds, esym)...)
  end
  unique(vcat(inds...))
end

"""    boundary_edges(ds)

Compute the edges of a 1D simplicial set that are either incident to in-degree 1 or out-degree 1 nodes.
For a graph, these are boundary vertices meaning leaf nodes. For our pipeflow problems,
these are the edges where material can enter the pipe network.
"""
function boundary_edges(ds)
  out_degree(x) = length(incident(ds, x, :∂v1))
  in_degree(x) = length(incident(ds, x, :∂v0))
  bpoints = findall(x -> out_degree(x) == 0 || in_degree(x) == 0, 1:nv(ds))
  sedges = vcat(incident(ds,bpoints,:∂v0)...)
  tedges = vcat(incident(ds,bpoints,:∂v1)...)
  bedges = collect(union(sedges, tedges))
  return bedges
end

"""    mask_boundary_edges(ds)

Provides the `boundary_edges(ds)` as a vector of 0/1 entries to use as a mask.
"""
function mask_boundary_edges(ds)
  D = ones(Int, ne(ds))
  D[boundary_edges(ds)] .= 0
  return D
end