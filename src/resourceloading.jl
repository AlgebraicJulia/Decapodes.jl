# 

struct Resource{T}
  path::Vector{String}
end 

abstract type ResourceProvider end

struct DirSsfsBackend <: ResourceProvider
  dir::String
end 

struct HDF5Backend <: ResourceProvider
  path::String
end 

# csv and stuff

load(p::ResourceProvider, r::Resource) = 

function load(p::DirSsfsBackend, r::Resource{Vector{Float64}})
  # TODO safe!
  path = p.dir * join(r.path, "/")
  content = read(path, String)
  map(s -> parse(Float64, s), split(content, " "))
end

function save(p::DirSsfsBackend, r::Resource{Vector{Float64}}, x::Vector{Float64})
  file = join(x, " ")
  write(p.path * join(r.path, "/"), file)
end
