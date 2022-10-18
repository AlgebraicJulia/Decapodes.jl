using Artifacts
using CombinatorialSpaces
using FileIO

#export AbstractMeshKey, loadmesh, UnitIcosphere, ThermoIcosphere, UnitUVSphere, ThermoUVSphere

abstract type AbstractMeshKey end

struct Icosphere{N, R} <: AbstractMeshKey
  n::N
  r::R
end

Icosphere(n) = Icosphere(n, 1.0)

function loadmesh(s::Icosphere)
  1 <= s.n <= 5 || error("The only icosphere divisions supported are 1-5")
  m = loadmesh_helper("UnitIcosphere$(s.n).obj")
  m[:point] .*= s.r
  return m
end

#loadmesh(meshkey::AbstractMeshKey)::EmbeddedDeltaSet2D

loadmesh_helper(obj_file_name) = EmbeddedDeltaSet2D(
  joinpath(artifact"all_meshes2", obj_file_name))

loadmesh(s, subdivision=Circumcenter()) = subdivide_duals!(loadmesh(s), subdivision)