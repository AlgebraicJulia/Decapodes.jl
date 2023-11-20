module Decapodes

using Catlab
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Programs
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.ACSetInterface
using MLStyle
using Base.Iterators

# TODO: Matt - Delete
import Unicode

using DiagrammaticEquations

export
# DerivOp,
# append_dot,
# normalize_unicode, 
findname, flat_op,
gensim, evalsim, closest_point, findnode, compile, compile_env,  PhysicsState,  default_dec_matrix_generate, default_dec_generate, default_dec_generate, VectorForm,
## meshes
AbstractMeshKey, loadmesh, Icosphere, Rectangle_30x10, Torus_30x10, Point_Map,
## coordinates 
CartesianPoint, SpherePoint, r, theta, phi, TangentBasis, θhat, ϕhat  

# TODO - Matt: Delete
# normalize_unicode(s::String) = Unicode.normalize(s, compose=true, stable=true, chartransform=Unicode.julia_chartransform)
# normalize_unicode(s::Symbol)  = Symbol(normalize_unicode(String(s)))
# DerivOp = Symbol("∂ₜ")
# TODO - Matt: this stays, but DE still has it
append_dot(s::Symbol) = Symbol(string(s)*'\U0307')

# include("migrate/decapodeacset.jl")
# include("migrate/language.jl")
# include("migrate/composition.jl")
# include("migrate/collages.jl")
include("coordinates.jl")
# include("migrate/visualization.jl")
include("simulation.jl")
include("meshes.jl")
# include("migrate/rewrite.jl")

end
