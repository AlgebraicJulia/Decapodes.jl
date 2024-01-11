using HDF5
using Mat

struct Resource{T}
  path::String
end 

abstract type ResourceProvider end


## SPACE-DELIMITED FILES
struct DirSsfsBackend <: ResourceProvider
  dir::String
end 

function load(p::DirSsfsBackend, r::Resource{Vector{Float64}})
  path = joinpath(p.dir, r.path) 
  content = read(path, String)
  map(s -> parse(Float64, s), split(content, " "))
end

function save(p::DirSsfsBackend, r::Resource{Vector{Float64}}, x::Vector{Float64})
  file = join(x, " ")
  write(joinpath(p.dir, r.path), file)
end

## CHARACTER-DELIMITED FILES
struct DirCdfsBackend <: ResourceProvider
  dir::String
  char::String
end

function load(p::DirCdfsBackend, r::Resource{Vector{Float64}})
  path = joinpath(p.dir, r.path) 
  content = h5open(path, "cw")
  map(s -> parse(Float64, s), split(content, r.char))
end

function save(p::DirCdfsBackend, r::Resource{Vector{Float64}}, x::Vector{Float64})
  file = join(x, r.char)
  write(joinpath(p.dir, r.path), file)
end


## HDF5 FILES
struct HDF5Backend <: ResourceProvider
  path::String
end 
## MAT
struct MatBackend <: ResourceProvider
  path::String
end

# load(p::ResourceProvider, r::Resource) = 




function load(p::HDF5Backend, r::Resource{Vector{Float64}})
  path = joinpath(p.dir, r.path)
  content = read(path, String)

end
#

## MAT
function load(p::MatBackend, r::Resource{Vector{Float64}})
  path = joinpath(p.dir, r.path)
  file = MAT.matread(path)
end

p=DirSsfsBackend("/home/you/Documents/UFAJ/AlgebraicJulia/Decapodes.jl/src/")
r=Resource{Vector{Float64}}("testssf.txt")
