using Artifacts
using CombinatorialSpaces
using FileIO

#export AbstractMeshKey, loadmesh, UnitIcosphere, ThermoIcosphere, UnitUVSphere, ThermoUVSphere

abstract type AbstractMeshKey end

#struct UnitIcosphere   <: AbstractMeshKey end
#struct ThermoIcosphere <: AbstractMeshKey end
#struct UnitUVSphere    <: AbstractMeshKey end
#struct ThermoUVSphere  <: AbstractMeshKey end

struct UnitIcosphere1 <: AbstractMeshKey end
struct UnitIcosphere2 <: AbstractMeshKey end
struct UnitIcosphere3 <: AbstractMeshKey end
struct UnitIcosphere4 <: AbstractMeshKey end
struct UnitIcosphere5 <: AbstractMeshKey end


#loadmesh(meshkey::AbstractMeshKey)::EmbeddedDeltaSet2D

loadmesh_helper(obj_file_name) = EmbeddedDeltaSet2D(
  joinpath(artifact"all_meshes2", obj_file_name))

#function loadmesh(meshkey::UnitIcosphere)::EmbeddedDeltaSet2D
#  all_meshes = artifact"all_meshes"
#  this_mesh = "unit_icosphere.obj"
#  mesh_obj = FileIO.load(File{format"OBJ"}(joinpath(all_meshes, this_mesh)))
#  mesh = EmbeddedDeltaSet2D(mesh_obj)
#end

function loadmesh(meshkey::UnitIcosphere1)::EmbeddedDeltaSet2D
  loadmesh_helper("UnitIcosphere1.obj")
end

function loadmesh(meshkey::UnitIcosphere2)::EmbeddedDeltaSet2D
  loadmesh_helper("UnitIcosphere2.obj")
end

function loadmesh(meshkey::UnitIcosphere3)::EmbeddedDeltaSet2D
  loadmesh_helper("UnitIcosphere3.obj")
end

function loadmesh(meshkey::UnitIcosphere4)::EmbeddedDeltaSet2D
  loadmesh_helper("UnitIcosphere4.obj")
end

function loadmesh(meshkey::UnitIcosphere5)::EmbeddedDeltaSet2D
  loadmesh_helper("UnitIcosphere5.obj")
end

