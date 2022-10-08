using Artifacts
using CombinatorialSpaces
using FileIO

#export AbstractMeshKey, loadmesh, UnitIcosphere, ThermoIcosphere, UnitUVSphere, ThermoUVSphere

abstract type AbstractMeshKey end

struct UnitIcosphere   <: AbstractMeshKey end
struct ThermoIcosphere <: AbstractMeshKey end
struct UnitUVSphere    <: AbstractMeshKey end
struct ThermoUVSphere  <: AbstractMeshKey end

#loadmesh(meshkey::AbstractMeshKey)::EmbeddedDeltaSet2D

loadmesh_helper(obj_file_name) = EmbeddedDeltaSet2D(
  joinpath(artifact"all_meshes", obj_file_name))

#function loadmesh(meshkey::UnitIcosphere)::EmbeddedDeltaSet2D
#  all_meshes = artifact"all_meshes"
#  this_mesh = "unit_icosphere.obj"
#  mesh_obj = FileIO.load(File{format"OBJ"}(joinpath(all_meshes, this_mesh)))
#  mesh = EmbeddedDeltaSet2D(mesh_obj)
#end
function loadmesh(meshkey::UnitIcosphere)::EmbeddedDeltaSet2D
  loadmesh_helper("unit_icosphere.obj")
end

function loadmesh(meshkey::ThermoIcosphere)::EmbeddedDeltaSet2D
  loadmesh_helper("thermo_icosphere.obj")
end

function loadmesh(meshkey::UnitUVSphere)::EmbeddedDeltaSet2D
  loadmesh_helper("unit_tie_gcm.obj")
end

function loadmesh(meshkey::ThermoUVSphere)::EmbeddedDeltaSet2D
  loadmesh_helper("thermo_tie_gcm.obj")
end

