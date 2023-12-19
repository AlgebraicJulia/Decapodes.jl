module Decapodes

using Catlab
using Catlab.Theories
# import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Programs
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.ACSetInterface
using MLStyle
using Base.Iterators

using DiagrammaticEquations
using DiagrammaticEquations.Deca

export
findname, flat_op,
gensim, evalsim, closest_point, findnode, compile, compile_env,  PhysicsState,  default_dec_matrix_generate, default_dec_generate, default_dec_generate, VectorForm,
AbstractMeshKey, loadmesh, Icosphere, Rectangle_30x10, Torus_30x10, Point_Map,
CartesianPoint, SpherePoint, r, theta, phi, TangentBasis, θhat, ϕhat  

append_dot(s::Symbol) = Symbol(string(s)*'\U0307')

include("coordinates.jl")
include("simulation.jl")
include("meshes.jl")

end
