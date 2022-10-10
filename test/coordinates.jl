using GeometryBasics
using Decapodes
using LinearAlgebra

x = CartesianPoint(Point3(1.0,2,3))
phi(x)
theta(x)
SpherePoint(x)
TangentBasis(x)
TangentBasis(x)(Point3(1.0,-1.0, 0)) == TangentBasis(x)(1.0,-1.0, 0)
TangentBasis(SpherePoint(x))(1,2,3) == Point3(1,2,3)