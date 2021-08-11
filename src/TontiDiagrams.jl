""" Tonti Diagrams

This file includes a framework for constructing and evaluating Tonti diagrams.
Tonti diagrams are stored as ACSets, and have an imperative interface for
describing physical variables and the relationships between them. This tooling
also lets a Tonti diagram be converted to a vectorfield, allowing for
simulation of physical systems through any traditional vectorfield solver.
"""

module TontiDiagrams

using Catlab.Present
using Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Graphs
import Catlab.Graphs: Graph
using Catlab.Graphics
using Catlab.Programs
using Catlab.WiringDiagrams
using CombinatorialSpaces
using CombinatorialSpaces: ⋆
using DifferentialEquations

export TheoryTontiDiagram, TontiDiagram, Space,
        add_variable!, add_variables!,
        add_derivative!, add_derivatives!,
        add_time_dep!,
        add_laplacian!, add_transition!, add_bc!,
        vectorfield, Open, gen_form


""" ACSet definition for a Tonti diagram.

This diagram is visualized in the q.uiver framework
[here](https://q.uiver.app/?q=WzAsMTAsWzEsMiwiTyJdLFsyLDIsIkkiXSxbMSwxLCJWIl0sWzEsMywiVCJdLFswLDIsIkJDIl0sWzAsMCwiVEQiXSxbMCwzLCJGdW5jIl0sWzMsMCwiQ29tcCJdLFszLDEsIkRpbWVuc2lvbiJdLFszLDIsIkxhYmVsIl0sWzQsMiwiYmN2Il0sWzUsMiwiaW50ZWciLDAseyJvZmZzZXQiOi0xfV0sWzAsMiwib3YiLDFdLFswLDMsIm90IiwxXSxbMSwyLCJpdiIsMV0sWzEsMywiaXQiLDFdLFs1LDIsImRlcml2IiwyLHsib2Zmc2V0IjoxfV0sWzQsNiwiYmNmdW5jIiwyLHsic3R5bGUiOnsiYm9keSI6eyJuYW1lIjoiZGFzaGVkIn19fV0sWzMsNiwidGZ1bmMiLDAseyJzdHlsZSI6eyJib2R5Ijp7Im5hbWUiOiJkYXNoZWQifX19XSxbMiw5LCJzeW1ib2wiLDEseyJzdHlsZSI6eyJib2R5Ijp7Im5hbWUiOiJkYXNoZWQifX19XSxbMiw4LCJkaW1lbnNpb24iLDEseyJzdHlsZSI6eyJib2R5Ijp7Im5hbWUiOiJkYXNoZWQifX19XSxbMiw3LCJjb21wbGV4IiwxLHsic3R5bGUiOnsiYm9keSI6eyJuYW1lIjoiZGFzaGVkIn19fV1d).
This Tonti diagram definition is currently very close to the simulation level
and generalized away from the DEC tooling. Future work will be done to provide
a closer integration between Tonti diagrams and DEC.

See Catlab.jl documentation for a description of the @present syntax.
"""
@present TheoryTontiDiagram(FreeSchema) begin
  Func::Data
  Comp::Data
  Dimension::Data
  Label::Data

  V::Ob
  I::Ob
  O::Ob
  T::Ob
  BC::Ob
  TD::Ob

  iv::Hom(I,V)
  it::Hom(I,T)
  ov::Hom(O,V)
  ot::Hom(O,T)
  bcv::Hom(BC,V)
  deriv::Hom(TD,V)
  integ::Hom(TD,V)

  tfunc::Attr(T,Func)
  bcfunc::Attr(BC,Func)

  complex::Attr(V,Comp)
  dimension::Attr(V,Dimension)
  symbol::Attr(V, Label)
  tlabel::Attr(T, Label)
  bclabel::Attr(BC, Label)
end

""" Space

Structure for storing the computed values for the primal/dual complex, the
boundary operators, the hodge stars, the laplacian operators, and the lie
operators. This is the operator which provides the DEC connection for the Tonti
diagra tooling.

Accessing the boundary, hodge, and laplacian operators is slightly complicated
by using arrays in a 1-indexed system. The hodge and boundary operators are
accessed with two indices, where the first determines the complex (1 for primal
and 2 for dual) and the second determines the dimension (1 for 0-forms, 2 for
1-forms, 3 for 2-forms, etc.). Thus, the boundary operator from dual 1-forms to
dual 2-forms is given as:

```julia
sp = Space(s)
sp.boundary[2,2]
```

The other goal of this structure is to cache the computed operators, yet
optimally this caching should be perfomed during the `vectorfield` operation.
"""
struct Space
  s::EmbeddedDeltaSet2D
  sd::EmbeddedDeltaDualComplex2D
  boundary
  hodge
  laplacian
  lie
end

""" Space(s::EmbeddedDeltaSet3D)

Caclulates all of the values stored in the `Space` object for a given complex
`s`.
"""
function Space(s::EmbeddedDeltaSet2D{O, P}) where {O, P}
  sd = EmbeddedDeltaDualComplex2D{O, eltype(P), P}(s)
  subdivide_duals!(sd, Barycenter())

  boundary = Array{typeof(d(0,s)), 2}(undef, 2,2)
  boundary[[1,2,3,4]] .= [d(0,s), dual_derivative(0,sd), d(1,s), dual_derivative(1,sd)]

  hodge = Array{typeof(⋆(0,sd)), 2}(undef, 2,3)
  hodge[[1,2,3,4,5,6]] .= [⋆(0,sd), inv(⋆(2,sd)), ⋆(1,sd), -1 .* inv(⋆(1,sd)), ⋆(2,sd), inv(⋆(0,sd))]

  # TODO:
  # Add Laplacian and Lie operators here
  laplacian  = Array{typeof(hodge[1,2]*boundary[1,1]),1}(undef,3)
  laplacian[1] = hodge[2,3]*boundary[2,2]*hodge[1,2]*boundary[1,1]

  laplacian[2] = hodge[2,2]*boundary[2,1]*hodge[1,3]*boundary[1,2] .+
                 boundary[1,1]*hodge[2,3]*boundary[2,2]*hodge[1,2]
  laplacian[3] = boundary[1,2]*hodge[2,2]*boundary[2,1]*hodge[1,3]

  lie_op = [hodge[2,3]*boundary[2,2]*hodge[1,2], hodge[2,2]*boundary[2,1]*hodge[1,3]]
  lie = [(∂s,s,v)->(∂s .= boundary[1,2]*∧(Tuple{0,1},sd,s,bound[2,1]*v)),
         (∂u,u,v)->(∂u .= (star_1*∧(Tuple{1,0},sd,v,inv_star_0*cobound_1_2*u) .+ cobound_0_1*star_2*∧(Tuple{1,1},sd,v,inv_star_1*u))),

        ]
  Space(s, sd, boundary, hodge, laplacian, nothing)
end

# TODO:
# Make more petri-style simulation ACSet
#
# Make more Tonti-aware ACSet w/ support of multiple complexes

# These functions are necessary for defining a TontiDiagram data structure from
# the ACSet
const AbstractTontiDiagram = AbstractACSetType(TheoryTontiDiagram)
const TontiDiagram = ACSetType(TheoryTontiDiagram,index=[:iv, :it, :ov, :ot, :bcv],
                                                 unique_index=[:symbol])

# Define an interface for an OpenTontiDiagram which allows Tonti diagrams to
# be composed over their corners
const OpenTontiDiagramOb, OpenTontiDiagram = OpenACSetTypes(TontiDiagram,
                                                            :V)

""" Open(td::TontiDiagram, states::Symbol)

Generates an OpenTontiDiagram with cospan legs on variables defined by the
symbols included in `states`. This OpenTontiDiagram can then be composed with
other OpenTontiDiagrams over a pattern given by an undirected wiring diagram.

```julia
OpenTontiDiagram(td, :x, :v)
```
"""
Open(td, states...) = OpenTontiDiagram{Function, Bool, Int64, Symbol}(td,
                        map(v->FinFunction([incident(td, v, :symbol)], nparts(td,:V)), states)...)

""" TontiDiagram()

Initialize an empty TontiDiagram object.
"""
TontiDiagram() = TontiDiagram{Function, Bool, Int64, Symbol}()

""" vectorfield(td::AbstractTontiDiagram, sp::Space)

Generates a Julia function which calculates the vectorfield of the Tonti
diagram. The state of the system is defined by a single vector which is a
flattening of all state variables of the system. Thus, this function returns
both the indices of each variable in the state-vector along with the
vectorfield function itself.

The resulting function has a signature of the form `f!(du, u, p, t)` and can be
passed to the DifferentialEquations.jl solver package.
"""
function vectorfield(td, sp::Space)
  fv_mem, ft_mem = init_mem(td, sp)
  dg = Graph(td)

  order = topological_sort(dg)

  state_vars = filter(x->length(incident(td,x,:ov))==0, 1:nparts(td, :V))

  input_vars = Dict{Symbol, Tuple{Int64,Int64}}()
  cur_head = 1
  for i in state_vars
    v_size = length(fv_mem[i])
    input_vars[td[i,:symbol]] = (cur_head,cur_head+v_size-1)
    cur_head += v_size
  end

  data_buffer = Dict{Type, Tuple{Array, Array}}(Float64=>(fv_mem, ft_mem))
  function system(du, u, t, p)
    if !(eltype(u) ∈ keys(data_buffer))
      data_buffer[eltype(u)] = init_mem(td, sp, type=eltype(u))
    end

    v_mem, t_mem = data_buffer[eltype(u)]

    for cur in order
      # Check if current is a variable or transition
      if cur > nparts(td, :V)
        cur -= nparts(td, :V)
        inputs = td[incident(td, cur, :it), :iv]
        outputs = incident(td, cur, :ot)
        td[cur,:tfunc](t_mem[outputs]..., v_mem[inputs]...)
      else
        inputs = incident(td, cur, :ov)
        if(length(inputs) == 0)
          # This means this is a state variable
          data_source = input_vars[td[cur,:symbol]]
          v_mem[cur] .= u[data_source[1]:data_source[2]]
        else
          v_mem[cur] .= 0
          for i in inputs
            v_mem[cur] .+= t_mem[i]
          end
        end
        bcs = incident(td, cur, :bcv)
        for bc in bcs
          td[bc, :bcfunc](v_mem[cur])
        end
      end
    end
    # If a state variable does not have a derivative defined, we keep it out
    # (we'll want to move these to the parameter argument instead)
    du .= 0
    for i in state_vars
      state_range = input_vars[td[i,:symbol]]
      out_var = td[incident(td, i, :integ), :deriv]
      if length(out_var) != 0
        du[state_range[1]:state_range[2]] .= v_mem[out_var[1]]
      end
    end
  end
  input_vars, system
end

""" add_variable!(td:TontiDiagram, symbol::Symbol, dimension::Int64, complex::Bool)

Adds a variable to the TontiDiagram system which can later be referenced by its
`symbol`. This constructor requires the dimensionality of the variable (0 ->
point, 1 -> line, etc.) and the complex it is defined on
(true -> primal/straight, false -> dual/twisted).

Definiting a system with a variable `v` defined on the primal lines would be
constructed as:
```@example
td = TontiDiagram()
add_variable!(td, :v, 1, true)
```
"""
function add_variable!(td, symbol::Symbol, dimension::Int64, complex::Bool)
  add_part!(td, :V,symbol=symbol, dimension=dimension, complex=complex)
end

""" add_variables!(td:TontiDiagram, vars::Tuple{Symbol, Int64, Bool}...)

Adds multiple variables to the TontiDiagram system which can later be
referenced by their symbols. This constructor follows the same pattern as
`add_variable!` with each variable specified as a tuple of:
```julia
(symbol, dimension, complex)
```

Defining a system with a variable `v` defined on the primal lines, a variable
`p̃` defined on the dual surfaces, and a variable `C` defined on the primal
surfaces would be constructed as:
```@example
td = TontiDiagram()
add_variables!(td, (:v, 1, true), (:p̃, 1, true), (:C, 2, true))
```
"""
function add_variables!(td, vars::Tuple{Symbol, Int64, Bool}...)
  for v in vars
    add_variable!(td, v...)
  end
end

""" add_transition!(td, dom_sym::Vector{Symbol}, func!, codom_sym::Vector{Symbol})

Adds a transition function from variables `dom_sym` to variables `codom_sym`
with its transition defined by `func!`. `func!` is expected to have the signature
`func(codom_sym..., dom_sym...)` and is expected to modify the values of the
`codom_sym` variables.

Defining a transition from variables `x` and `y` to `m` that calculates the
magnitude of the values in `x` and `y` as vector coordinates would be defined as:

```@example
td = TontiDiagram()
add_variables!(td, (:x, 0, 1), (:y, 0, 1). (:m, 0, 1))
add_transition!(td, [:x, :y], (m,x,y)->(m .= sqrt.(x .* y)), [:m])
```
"""
function add_transition!(td, dom_sym::Array{Symbol,1}, func!, codom_sym::Array{Symbol,1}; label=Symbol(""))
  dom = [findfirst(v->v == s, td[:symbol]) for s in dom_sym]
  codom = [findfirst(v->v == s, td[:symbol]) for s in codom_sym]
  t = add_part!(td, :T, tfunc=func!, tlabel=label)
  add_parts!(td, :I, length(dom), iv=dom, it=t)
  add_parts!(td, :O, length(codom), ov=codom, ot=t)
end

""" add_derivative!(td::TontiDiagram, sp::Space, dom_sym::Symbol,
codom_sym::Symbol)

Adds a derivative transition from variable `dom_sym` to variable `codom_sym`
using the boundary operators from `sp`. This function determines which boundary
operator to use and inserts an appropriate transition between the two
variables.

Defining a spatial derivative relationship between the primal 0-form `x` and
the primal 1-form `Δx` cis given as follows:
```julia
add_derivative(td, sp, :x, :Δx)
```
"""
function add_derivative!(td, sp, dom_sym, codom_sym)
  dom = findfirst(v->v==dom_sym, td[:symbol])
  codom = findfirst(v->v==codom_sym, td[:symbol])

  # TODO:
  # Add tests for proper dimensions, complexes, etc.
  # This will later be replaced as we pre-initialize all boundary operators
  bound = sp.boundary[(td[dom,:complex] ? 1 : 2), td[dom,:dimension]+1]
  func(x,y) = (x.=bound*y)
  add_transition!(td, [dom_sym],func,[codom_sym]; label=:d)
end

""" add_derivatives!(td::TontiDiagram, sp::Space, vars:Pair{Symbol, Symbol}...)

Adds multiple derivative transition between pairs of variables, using the same
syntax as in `add_derivative!`.

Example usage:
```julia
add_derivatives!(td, sp, (:x,:Δx), (:y, :Δy))
```
"""
function add_derivatives!(td, sp, vars::Pair{Symbol, Symbol}...)
  for v in vars
    add_derivative!(td, sp, v[1],v[2])
  end
end

""" add_time_dep!(td::TontiDiagram, deriv_sym::Symbol, integ_sym::Symbol)

Adds a time derivative relationship between the variables `deriv_sym` and
`integ_sym` (where `deriv_sym` is the time derivative of `integ_sym`). These
relationships are used to determine the state-variables of the system.
"""
function add_time_dep!(td, deriv_sym::Symbol, integ_sym::Symbol)
  deriv = findfirst(v->v==deriv_sym, td[:symbol])
  integ = findfirst(v->v==integ_sym, td[:symbol])

  add_part!(td, :TD, integ=integ, deriv=deriv)
end

""" add_bc!(td::TontiDiagram, var_sys::Symbol, func::Function)

Adds a "boundary condition" to the variable `var_sys` by applying `func` to the
values of this variable during simulation. This function is the last one
evaluated on the data of the variable `var_sys`, and so can be used to enforce
any relevant boundary conditions.

TODO:
Add time dependency of boundary condition function to allow for time-varying BCs.
"""
function add_bc!(td, var_sym, func; label=:BC)
  var = findfirst(v->v==var_sym, td[:symbol])
  add_part!(td, :BC, bcfunc=func, bcv=var, bclabel=label)
end

# Note: This function can be made more efficient if combined with existing
# transformations.
# e.g. Advection-diffusion can be merged after the initial wedge
# product/coboundary operator
#
# Currently only defined on primal complices (can this be applied to dual
# complices?)

""" add_laplacian!(td::TontiDiagrams, sp::Space, dom_sym::Symbol, codom_sym::Symbol; coef::Float64)

Adds a transition which defines `codom_sym` as the laplacian of `dom_sym` with
a constant scaling factor of `coef`.
"""
function add_laplacian!(td, sp, dom_sym, codom_sym; coef=1.0)
  sd = sp.sd
  dom = findfirst(v->v==dom_sym, td[:symbol])
  codom = findfirst(v->v==codom_sym, td[:symbol])

  lap_op = sp.laplacian[td[dom,:dimension]+1] # laplace_beltrami(Val{td[dom,:dimension]},sd)
  func(x,y) = (x .= coef * (lap_op*y))
  add_transition!(td, [dom_sym], func, [codom_sym]; label=:Δ)
end

function init_mem(td, s::EmbeddedDeltaSet1D)
  # Fill out this function
end

function init_mem(td, sp::Space; type=Float64)
  s = sp.s
  primal_size = [nv(s), ne(s), ntriangles(s)]
  dual_size   = [ntriangles(s), ne(s), nv(s)]

  t_mem = Array{Array{type,1},1}()
  v_mem = Array{Array{type,1},1}()

  for i in 1:nparts(td, :O)
    var = td[i,:ov]
    push!(t_mem, zeros(type,
                      td[var,:complex] ? primal_size[td[var,:dimension]+1] :
                      dual_size[td[var,:dimension]+1]))
  end

  for v in 1:nparts(td,:V)
    push!(v_mem, zeros(type,
                      td[v,:complex] ? primal_size[td[v,:dimension]+1] :
                      dual_size[td[v,:dimension]+1]))
  end

  v_mem, t_mem
end

function Graph(td; extended=false)
  g = Graph()
  add_vertices!(g, nparts(td, :V) + nparts(td, :T))
  nvars = nparts(td, :V)
  for i in 1:nparts(td, :I)
    add_edge!(g, td[i,:iv], td[i,:it] + nvars)
  end
  for o in 1:nparts(td, :O)
    add_edge!(g, td[o,:ot] + nvars, td[o,:ov])
  end
  if extended
    for bc in 1:nparts(td, :BC)
      add_edge!(g, td[bc, :bcv], td[bc, :bcv])
    end
    for t in 1:nparts(td, :TD)
      add_edge!(g, td[t, :deriv], td[t, :integ])
    end
  end
  g
end

""" Construct a 0-form based on a scalar function

This operator accepts a scalar function and evaulates it at each point on the
simplex, returning a 0-form.
"""
function gen_form(s::EmbeddedDeltaSet2D, f::Function)
  map(f, point(s))
end


function propertygraph(td::TontiDiagram;
                       prog::AbstractString="neato", graph_attrs::AbstractDict=Dict(),
                       node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(:len=>"1.5"),
                       node_labels::Bool=true, edge_labels::Bool=false)

  nvars = nparts(td, :V)
  isVar(v) = (v <= nvars)
  g = Graph(td, extended = true)

  dim_colors = ["green", "yellow", "red", "blue"]

	node_labeler(v) = begin
    if isVar(v)
      Dict(:label=>"$(td[v, :symbol])",
           :style=>td[v,:complex] ? "filled" : "filled, dotted",
           :width=>"0.75", :height=>"0.75", :fixedsize=>"false",
           :shape=>"circle",
           :fillcolor=>dim_colors[td[v,:dimension]+1],
           )
    else
      Dict(:label=>"$(td[v - nvars, :tlabel])",
           :width=>"0.75", :height=>"0.75", :fixedsize=>"false",
           :shape=>"square"
           )
    end
  end

  edge_labeler(e) = begin
    if e <= nparts(td, :I) + nparts(td, :O)
      Dict(:color=>"black")
    else
      e -= nparts(td, :I) + nparts(td, :O)
      if e <= nparts(td, :BC)
        Dict(:color=>"black",
             :label=>"$(td[e, :bclabel])")
      else
        e -= nparts(td, :BC)
        Dict(:color=>"black", :style=>"dotted")
      end
    end
  end


  PropertyGraph{Any}(g, node_labeler, edge_labeler;
                        prog = prog,
                        graph = merge!(Dict(:rankdir => "TB"), graph_attrs),
                        node = merge!(Graphics.GraphvizGraphs.default_node_attrs(node_labels), node_attrs),
                        edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
                        )
end

function Graphics.to_graphviz(td::TontiDiagram; kwargs...)
  propertygraph(td; kwargs...) |> to_graphviz
end

end
