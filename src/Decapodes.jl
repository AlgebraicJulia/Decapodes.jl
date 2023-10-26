module Decapodes

using Catlab
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Programs
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using Catlab.ACSetInterface
using LinearAlgebra
using MLStyle
using Base.Iterators
using SparseArrays

import Unicode

export normalize_unicode, DerivOp, append_dot, unicode!, vec_to_dec!,
  SchDecapode, SchNamedDecapode, AbstractDecapode, AbstractNamedDecapode, Decapode, NamedDecapode, SummationDecapode, fill_names!, dot_rename!, expand_operators, add_constant!, add_parameter, infer_types!, resolve_overloads!, op1_inf_rules_1D, op2_inf_rules_1D, op1_inf_rules_2D, op2_inf_rules_2D, op1_res_rules_1D, op2_res_rules_1D, op1_res_rules_2D, op2_res_rules_2D,
  Term, Var, Judgement, Eq, AppCirc1, AppCirc2, App1, App2, Plus, Tan, term, parse_decapode,
  VectorForm, PhysicsState, findname, findnode,
  compile, compile_env, gensim, evalsim, closest_point, flat_op,
  AbstractMeshKey, loadmesh, Icosphere, Rectangle_30x10, Torus_30x10, Point_Map,
  Open, OpenSummationDecapodeOb, OpenSummationDecapode, unique_by, unique_by!, oapply,
  CartesianPoint, SpherePoint, r, theta, phi, TangentBasis, θhat, ϕhat,
  average_rewrite, recursive_delete_parents, contract_operators,
  default_dec_matrix_generate, default_dec_generate, default_dec_generate_1D, default_dec_generate_2D,
  @decapode

normalize_unicode(s::String) = Unicode.normalize(s, compose=true, stable=true, chartransform=Unicode.julia_chartransform)
normalize_unicode(s::Symbol)  = Symbol(normalize_unicode(String(s)))
DerivOp = Symbol("∂ₜ")
append_dot(s::Symbol) = Symbol(string(s)*'\U0307')

include("decapodeacset.jl")
include("language.jl")
include("composition.jl")
include("coordinates.jl")
include("visualization.jl")
include("operators.jl")
include("simulation.jl")
include("meshes.jl")
include("rewrite.jl")

end
