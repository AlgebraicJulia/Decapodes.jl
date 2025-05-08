using CombinatorialSpaces
using Catlab
using DiagrammaticEquations
using Decapodes
using Decapodes.Canon.Physics: oscillator
using GeometryBasics: Point3
Point3D = Point3{Float64}

mesh1 = EmbeddedDeltaSet1D{Bool, Point3D}()
mesh2 = EmbeddedDeltaSet1D{Bool, Point3D}()
dualmesh1 = EmbeddedDeltaDualComplex1D{Bool, Float64, Point3D}(mesh1)
dualmesh2 = EmbeddedDeltaDualComplex1D{Bool, Float64, Point3D}(mesh2)
subdivide_duals!(dualmesh1, Circumcenter())
subdivide_duals!(dualmesh2, Circumcenter())

D = FreeDiagram(FinSet.([3,2,3]), # list of objects
 [ # list of (hom, src, tgt) tuples
  (FinFunction([1,2], 3), 2,1),
  (FinFunction([1,2], 3), 2,3),
  ]
)

(→)(x::ACSet, y::ACSet) = ACSetTransformation(x, y)

J = FreeDiagram([dualmesh1, dualmesh2],
            [(dualmesh1 → dualmesh2, 1, 2),
             (dualmesh1 → dualmesh2, 1, 2)])
