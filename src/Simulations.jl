""" Generating simulation functions from directed wiring diagrams

This module currently provides support for generating a simple explicit
time-stepping function based on a directed wiring diagram. The generated
function is compatible with the DifferentialEquations.jl interface, which is
used in demo execution.
"""

module Simulations

using CombinatorialSpaces
using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra
using Catlab.Theories
using PreallocationTools
using LinearAlgebra

export gen_sim,
       BoxFunc, MatrixFunc, ElementwiseFunc, ArbitraryFunc, InPlaceFunc,
       ConstantFunc, TDInPlaceFunc

abstract type BoxFunc end

""" MatrixFunc

Linear operator applied through matrix multiplication.
"""
struct MatrixFunc <: BoxFunc end

""" ElementwiseFunc

Function applied to each element of one or more vectors.
"""
struct ElementwiseFunc <: BoxFunc end

""" ArbitraryFunc

Out of place function which applies to any number of arguments.
"""
struct ArbitraryFunc <: BoxFunc end

""" ConstantFunc

Function which always returns the same value.
"""
struct ConstantFunc <: BoxFunc end

""" InPlaceFunc

Function which operates on arguments in-place, not returning a value.
"""
struct InPlaceFunc <: BoxFunc end

""" TDInPlaceFunc

Function which accepts the current time as the last argument. Used to set up
time-dependent physics or boundary conditions.
"""
struct TDInPlaceFunc <: BoxFunc end

form2dim = Dict(:Scalar => x->1,
                :Form0 => nv,
                :Form1 => ne,
                :Form2 => ntriangles,
                :DualForm2 => nv,
                :DualForm1 => ne,
                :DualForm0 => ntriangles)
dims(x) = begin
  k = filter(i -> x[:type] isa ObExpr{i}, keys(form2dim))
  length(k) == 1 || error("Object $x has multiple keys in `form2dim`")
  form2dim[first(k)]
end

""" gen_sim(dwd::WiringDiagram, name2func::Dict{Symbol, <:BoxFunc},
            s::EmbeddedDeltaDualComplex2D;
            form2dim=form2dim, params=[], autodiff=false)

This function generates a function which evaluates the acylcic directed wiring
diagram `dwd`, using the functions in `name2func` for operators to execute
corresonding to each box and information from the mesh `s` to pre-allocate
necessary memory. This operator can generate a function which is compatible
with the autodifferentiation solvers in DifferentialEquations.jl.
"""
function gen_sim(dwd::WiringDiagram, name2func, s; form2dim=form2dim, params=[], autodiff=false)
  check_consistency(dwd)
  d = dwd.diagram

  # Generate cached memory. These variables are indexed by their respective
  # output ports
  mem = [zeros(Float64, dims(f)(s)) for f in d[:out_port_type]]
  if autodiff
    mem = [dualcache(zeros(Float64, dims(f)(s))) for f in d[:out_port_type]]
  end
  tgt2src = Dict{Int64, Int64}()
  for w in 1:nparts(d, :Wire)
    d[w, :tgt] ∈ keys(tgt2src) && error("Two wires input to port $(d[w, :src])")
    tgt2src[d[w, :tgt]] = d[w, :src]
  end

  # Get the indices in input/output arrays for specific input/output forms
  input_inds = Vector{Tuple{Int64, Int64, Bool}}()
  output_inds = Vector{Tuple{Int64, Int64}}()

  cur_ind = 0
  cur_param = 0
  for ft in d[:outer_in_port_type]
    mag = dims(ft)(s)
    if ft[:name] ∈ params
      push!(input_inds, (cur_param+1, mag + cur_param, true))
      cur_param += mag
    else
      push!(input_inds, (cur_ind+1, mag + cur_ind, false))
      cur_ind += mag
    end
  end
  cur_ind = 0
  for ft in d[:outer_out_port_type]
    mag = dims(ft)(s)
    push!(output_inds, (cur_ind+1, mag + cur_ind))
    cur_ind += mag
  end

  # Make local variables for any MatrixFunc/ArbitraryFunc functions to be
  # included in the function closure
  n2f = deepcopy(name2func)
  matrices = Vector{Any}()
  funcs = Vector{Function}()

  # First develop a list of needed functions/matrices
  used_funcs = unique(d[:value])

  for k in used_funcs
    ops = split("$k", "⋅")
    if length(ops) > 1
      mat_val = foldl(*, map(o->name2func[Symbol(o)][:operator], ops))
      push!(matrices, mat_val)
      n2f[k] = Dict(:operator => :(matrices[$(length(matrices))]), :type => MatrixFunc())
    else
      v = name2func[k]
      if(v[:type] isa MatrixFunc)
        push!(matrices, v[:operator])
        n2f[k] = Dict(:operator => :(matrices[$(length(matrices))]), :type => v[:type])
      elseif(v[:type] isa Union{ElementwiseFunc, ArbitraryFunc, InPlaceFunc, TDInPlaceFunc})
        push!(funcs, v[:operator])
        n2f[k] = Dict(:operator => :(funcs[$(length(funcs))]),:type => v[:type])
      end
    end
  end

  # Compute composed matrices
  for v in d[:value]

  end

  # Topological sort of DWD for scheduling
  execution_order = topological_sort(internal_graph(dwd))

  # Generate function expression
  body = quote
  end

  # Assign variables for each dual memory location
  # Allows for get_tmp to be called for each cached memory location
  if autodiff
    dual_init = map(1:length(mem)) do i
      :($(Symbol("d_$i")) = get_tmp(mem[$i], u))
    end
    append!(body.args, dual_init)
  end

  exec_dwd = map(execution_order) do b
    iw = in_wires(dwd, b)

    input_args = map(iw) do w
      p = w.source
      if p.box == -2
        if input_inds[p.port][3]
          :(view(p, $(input_inds[p.port][1]):$(input_inds[p.port][2])))
        else
          :(view(u, $(input_inds[p.port][1]):$(input_inds[p.port][2])))
        end
      else
        if autodiff
          :($(Symbol("d_$(incident(d, p.box, :out_port_box)[p.port])")))
        else
          :(mem[$(incident(d, p.box, :out_port_box)[p.port])])
        end
      end
    end
    input_args[[w.target.port for w in iw]] .= input_args

    output_args = map(incident(d, b, :out_port_box)) do p
      if autodiff
        :($(Symbol("d_$p")))
      else
        :(mem[$p])
      end
    end

    gen_func(input_args, output_args, n2f[d[b, :value]])
  end
  append!(body.args, exec_dwd)

  # set du values
  du_set = map(1:nparts(d, :OuterOutPort)) do i
    iw = incident(d, i, :out_tgt)
    length(iw) == 1 || error(Output)
    w = first(iw)
    pt = d[w, :out_tgt]
    ps = d[w, :out_src]
    if autodiff
      :(setindex!(du, Symbol("d_$ps"), $(collect(output_inds[pt][1]:output_inds[pt][2]))))
    else
      :(setindex!(du, mem[$ps], $(collect(output_inds[pt][1]:output_inds[pt][2]))))
    end
  end
  append!(body.args, du_set)

  # Currently this generates a function at the module (can be accessed from the module)
  # Is this alright, or is this going to just lead to memory problems?
  fname = gensym()
  res = :($fname(du, u, p, t, mem, matrices, funcs) = $body)
  eval(res)
  f_local = deepcopy(eval(fname))
  function sim_func(du, u, p, t)
    f_local(du, u, p, t, mem, matrices, funcs)
  end, res
end

function gen_func(input_args, output_args, func_info)
  _gen_func(func_info[:type], input_args, output_args, func_info[:operator])
end

function _gen_func(::MatrixFunc, input_args, output_args, func)
  length(input_args) == 1 || error("Length of input for Matrix function is $(length(input_args)) ≠ 1")
  length(output_args) == 1 || error("Length of output for Matrix function is $(length(output_args)) ≠ 1")
  :(mul!($(output_args[1]), $(func), $(input_args[1])))
end

function _gen_func(::ArbitraryFunc, input_args, output_args, func)

  length(output_args) == 0 && error("Length of output for Arbitrary function is $(length(output_args))")
  :($(Meta.parse(join(output_args, ","))) = $(Expr(:call, func, input_args...)))
end

function _gen_func(::InPlaceFunc, input_args, output_args, func)
  length(output_args) == 0 && error("Length of output for Arbitrary function is $(length(output_args))")
  :($(Expr(:call, func, vcat(output_args, input_args)...)))
end

function _gen_func(::TDInPlaceFunc, input_args, output_args, func)
  _gen_func(InPlaceFunc(), vcat(input_args, [:t]), output_args, func)
end

function _gen_func(::ConstantFunc, input_args, output_args, func)
  length(output_args) == 0 && error("Length of output for Constant function is $(length(output_args))")
  :($(Meta.parse(join(output_args, ","))) .= $(func))
end


function _gen_func(::ElementwiseFunc, input_args, output_args, func)

  length(output_args) == 0 && error("Length of output for Arbitrary function is $(length(output_args))")
  body = :($(Meta.parse(join([Expr(:ref, o, :i) for o in output_args], ","))) = $(Expr(:call, func, [Expr(:ref, n, :i) for n in input_args]...)))
  quote
    for i in eachindex($(input_args[1]))
      $body
    end
  end
end

function check_consistency(dwd::WiringDiagram)
  d = dwd.diagram
  for w in 1:nparts(d, :Wire)
    sp = d[w, :src]
    tp = d[w, :tgt]
    stype = d[sp, :out_port_type]
    ttype = d[tp, :in_port_type]
    sb = d[sp, :out_port_box]
    tb = d[tp, :in_port_box]

    stype == ttype ||
      error("""Wire $w maps a port of type $stype to a port of type $ttype
            This wire is between box $(sb)($(d[sb, :value])) and box $(tb)($(d[tb, :value]))""")
  end
end

end
