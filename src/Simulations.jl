module Simulations

using CombinatorialSpaces
using Catlab
using Catlab.WiringDiagrams
using Catlab.CategoricalAlgebra

export gen_sim,
       BoxFunc, MatrixFunc, ElementwiseFunc, ArbitraryFunc, ConstantFunc

abstract type BoxFunc end
struct MatrixFunc <: BoxFunc end
struct ElementwiseFunc <: BoxFunc end
struct ArbitraryFunc <: BoxFunc end
struct ConstantFunc <: BoxFunc end

form2dim = Dict(:Scalar => x->1,
        				:Form0 => nv,
        				:Form1 => ne,
        				:Form2 => ntriangles,
        				:DForm2 => nv,
        				:DForm1 => ne,
        				:DForm0 => ntriangles)

function gen_sim(dwd::WiringDiagram, name2func, s; form2dim=form2dim)
  check_consistency(dwd)
  d = dwd.diagram

  # Generate cached memory. These variables are indexed by their respective
  # output ports
  mem = [zeros(Float64, form2dim[f](s)) for f in d[:out_port_type]]
  tgt2src = Dict{Int64, Int64}()
  for w in 1:nparts(d, :Wire)
    d[w, :tgt] ∈ keys(tgt2src) && error("Two wires input to port $(d[w, :src])")
    tgt2src[d[w, :tgt]] = d[w, :src]
  end

  # Get the indices in input/output arrays for specific input/output forms
  input_inds = Vector{Tuple{Int64, Int64}}()
  output_inds = Vector{Tuple{Int64, Int64}}()

  cur_ind = 0
  for ft in d[:outer_in_port_type]
    mag = form2dim[ft](s)
    push!(input_inds, (cur_ind+1, mag + cur_ind))
    cur_ind += mag
  end
  cur_ind = 0
  for ft in d[:outer_out_port_type]
    mag = form2dim[ft](s)
    push!(output_inds, (cur_ind+1, mag + cur_ind))
    cur_ind += mag
  end

  # Make local variables for any MatrixFunc/ArbitraryFunc functions to be
  # included in the function closure
  n2f = deepcopy(name2func)
  matrices = Vector{Any}()
  funcs = Vector{Function}()
  for (k,v) in name2func
    if(v[:type] isa MatrixFunc)
      push!(matrices, v[:operator])
      n2f[k] = Dict(:operator => :(matrices[$(length(matrices))]), :type => v[:type])
    elseif(v[:type] isa Union{ElementwiseFunc, ArbitraryFunc})
      push!(funcs, v[:operator])
      n2f[k] = Dict(:operator => :(funcs[$(length(funcs))]),:type => v[:type])
    end
  end

  # Topological sort of DWD for scheduling
  execution_order = topological_sort(internal_graph(dwd))

  # Generate function expression
  body = quote
  end

  exec_dwd = map(execution_order) do b
    iw = in_wires(dwd, b)
    ow = out_wires(dwd, b)

    input_args = map(iw) do w
      p = w.source
      if p.box == -2
        :(u[$(input_inds[p.port][1]):$(input_inds[p.port][2])])
      else
        :(mem[$(incident(d, p.box, :out_port_box)[p.port])])
      end
    end
    input_args[[w.target.port for w in iw]] .= input_args

    output_args = map(ow) do w
      :(mem[$(incident(d, b, :out_port_box)[w.source.port])])
    end
    output_args[[w.source.port for w in ow]] .= output_args

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
    :(du[$(output_inds[pt][1]):$(output_inds[pt][2])] .= mem[$ps])
  end
  append!(body.args, du_set)

  res = :(f!(du, u, p, t, mem, matrices, funcs) = $body)
  eval(res)
  function sim_func(du, u, p, t)
    f!(du, u, p, t, mem, matrices, funcs)
  end, res
end

function gen_func(input_args, output_args, func_info)
  _gen_func(func_info[:type], input_args, output_args, func_info[:operator])
end

function _gen_func(::MatrixFunc, input_args, output_args, func)
  length(input_args) == 1 || error("Length of input for Matrix function is $(length(input_args)) ≠ 1")
  length(output_args) == 1 || error("Length of output for Matrix function is $(length(output_args)) ≠ 1")
  :($(output_args[1]) .= $(func) * $(input_args[1]))
end

function _gen_func(::ArbitraryFunc, input_args, output_args, func)

  length(output_args) == 0 && error("Length of output for Arbitrary function is $(length(output_args))")
  :($(Meta.parse(join(output_args, ","))) = $(Expr(:call, func, input_args...)))
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
