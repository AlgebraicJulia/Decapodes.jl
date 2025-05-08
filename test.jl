using Revise

using CombinatorialSpaces
using ACSets

s = loadmesh(Icosphere(8));

q = From(:Tri) |> Where([:∂e0, :∂e1, :∂e2], (x,y,z) -> x + y + z < 100) & Where(:Tri, iseven)
@time q(s)

out = Int64[]
function foo!(s::HasDeltaSet, out::AbstractVector{Int64})
    for tri in parts(s, :Tri)
        if iseven(tri)
            e0::Int64 = s[tri,:∂e0]
            e1::Int64 = s[tri,:∂e1]
            e2::Int64 = s[tri,:∂e2]
            if e0 + e1 + e2 < 100
                push!(out, tri)
            end
        end
    end
    out
end

@time foo!(s, Int64[])


