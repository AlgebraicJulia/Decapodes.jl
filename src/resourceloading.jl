# 

struct Resource
    path::Vector{String}
    type::Type
end 

abstract type ResourceProvider end

struct DirCsvsBackend <: ResourceProvider
    path::String
end 

struct HDF5Backend <: ResourceProvider
    path::String
end 

load(p::ResourceProvider, r::Resource)

save(p::ResourceProvider, r::Resource, x)
