using Catlab
using Catlab.Theories
import Catlab.Theories: otimes, oplus, compose, ⊗, ⊕, ⋅, associate, associate_unit, Ob, Hom, dom, codom
using Catlab.Present
using Catlab.Programs
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.DirectedWiringDiagrams
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using LinearAlgebra
using MLStyle
using Base.Iterators

import Unicode

export normalize_unicode, DerivOp, append_dot,
  SchDecapode, SchNamedDecapode, AbstractDecapode, AbstractNamedDecapode, Decapode, NamedDecapode, SummationDecapode, fill_names!, expand_operators,
  Term, Var, Judge, Eq, AppCirc1, AppCirc2, App1, App2, Plus, Tan, term, parse_decapode,
  VectorForm, PhysicsState, findname, findnode,
  compile, compile_env, gensim, closest_point, flat_op,
  Open, OpenSummationDecapodeOb, OpenSummationDecapode, unique_by, unique_by!, oapply

normalize_unicode(s::String) = Unicode.normalize(s, compose=true, stable=true, chartransform=Unicode.julia_chartransform)
normalize_unicode(s::Symbol)  = Symbol(normalize_unicode(String(s)))
DerivOp = Symbol("∂ₜ")
append_dot(s::Symbol) = Symbol(string(s)*'\U0307')

include("decapodeacset.jl")
include("language.jl")
include("composition.jl")
include("visualization.jl")
include("simulation.jl")
