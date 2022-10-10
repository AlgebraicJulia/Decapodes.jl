using GeometryBasics
using LinearAlgebra

"""    CartesianPoint{T}(p)

a point in cartesian coordinates, intended as a wrapper around Point3 from GeometryBasics.
"""
struct CartesianPoint{T}
    p::T
end

r(p::CartesianPoint) = hypot(p.p)
theta(p::CartesianPoint) = atan(sqrt(p.p[1]^2 + p.p[2]^2)/p.p[3])
phi(p::CartesianPoint) = atan(p.p[2]/p.p[1])
  
"""    SpherePoint{T}(p)

a point in spherical coordinates, intended as a wrapper around Point3 from GeometryBasics.
"""
struct SpherePoint{T}
    p::T
end

r(p::SpherePoint) = p.p[1]
theta(p::SpherePoint) = p.p[2]
phi(p::SpherePoint) = p.p[3]

SpherePoint(p::CartesianPoint{T}) where T = SpherePoint{T}(T(r(p), theta(p), phi(p)))


struct TangentBasis{P}
  p::P
end

θhat(p::TangentBasis{P} where {P<:SpherePoint}) = Point2(1.0, 0.0)
ϕhat(p::TangentBasis{P} where {P<:SpherePoint}) = Point2(0.0, 1.0)

θhat(b::TangentBasis{P} where {P<:CartesianPoint}) = Point3(cos(theta(b.p))*cos(phi(b.p)), cos(theta(b.p))*sin(phi(b.p)), -sin(theta(b.p)))*sign(b.p.p[3])
ϕhat(b::TangentBasis{P} where {P<:CartesianPoint}) = Point3(-sin(phi(b.p)), cos(phi(b.p)), 0)

"""    tb(w)

Take a linear combinations of the tangent vectors at the base point. 
Use this to get a vector tangent to the sphere in the coordinate system of the base point.
If the base point is in spherical coordinates, this is the identity,
if the base point is in cartesian coordinates, it returns the tangent vector in cartesian coordinates.
"""
function (tb::TangentBasis)(w) 
    tb(w[1], w[2])
end

function (tb::TangentBasis{P})(w) where {P <: SpherePoint}
    return w
end

function (tb::TangentBasis)(w1,w2) 
    w1p = w1*θhat(tb)
    w2p = w2*ϕhat(tb)
    return w1p + w2p
end