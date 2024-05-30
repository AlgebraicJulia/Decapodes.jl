using GeometryBasics
using Decapodes
using LinearAlgebra

x = CartesianPoint(Point3(1.0,2,3))
# TODO: Verify these are the correct literals.
@test isapprox(phi(x),   1.1071487, atol=1e-7)
@test isapprox(theta(x), 0.6405223, atol=1e-7)
@test isapprox(r(x),     3.7416573, atol=1e-7)

@test typeof(SpherePoint(x)) == SpherePoint{Point3{Float64}}
@test isapprox(SpherePoint(x).p[1], r(x), atol=1e-7)
@test isapprox(SpherePoint(x).p[2], theta(x), atol=1e-7)
@test isapprox(SpherePoint(x).p[3], phi(x), atol=1e-7)

@test typeof(TangentBasis(x)) == TangentBasis{CartesianPoint{Point3{Float64}}}
@test TangentBasis(x)(Point2(1.0,-1.0)) == TangentBasis(x)(1.0,-1.0)
@test TangentBasis(SpherePoint(x))(1,2) == Point2(1,2)