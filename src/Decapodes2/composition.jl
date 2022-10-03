
import Catlab.WiringDiagrams: oapply

"""      function unique_by!(acset, column_names::Vector{Symbol})

Given column names from the same table, remove duplicate rows.

WARNING: This function does not check if other tables index into the one given.
Removal of rows is performed with prejudice.

See also: [`unique_by`](@ref).

# Examples
```julia-repl
julia> unique_by!(parallel_arrows(Graph, 123), :E, [:src,:tgt]) == parallel_arrows(Graph, 1)
true
```
"""
function unique_by!(acset, table::Symbol, columns::Vector{Symbol})
  # Note: Declarative CT methods are prefered to imperative index arithmetic.
  rows = mapreduce(x -> acset[x], hcat, columns)
  rem_parts!(acset, table,
    setdiff(parts(acset, table),
      unique(i -> rows[i,:], eachindex(eachrow(rows)))))
  return acset
end

"""      function unique_by(acset, column_names::Vector{Symbol})

Given column names from the same table, return a copy of the acset with
duplicate rows removed. Removal of rows is performed with prejudice.

WARNING: This function does not check if other tables index into the one given.
Removal of rows is performed with prejudice.

See also: [`unique_by!`](@ref).

# Examples
```julia-repl
julia> unique_by(parallel_arrows(Graph, 123), :E, [:src,:tgt]) == parallel_arrows(Graph, 1)
true
```
"""
function unique_by(acset, table::Symbol, columns::Vector{Symbol})
  acset_copy = copy(acset)
  unique_by!(acset_copy, table, columns)
end

"""    function type_check_decapodes_composition(decapodes_vars, relation, local_ports)

Check that the types of all Vars connected by the same junction match.

This function only throws an error on the first type mismatch found.
"""
function type_check_decapodes_composition(decapodes_vars, relation, local_ports)
  r = relation
  decs = first.(decapodes_vars)
  vars = last.(decapodes_vars)
  for j ∈ junctions(r)
    # Get the type of the first variable attached to this junction.
    P = ports_with_junction(r, j)
    p₁ = first(P)
    b₁ = r[p₁, :box]
    lp₁ = local_ports[p₁]
    type₁ = decs[b₁][lp₁, :type]
    # Compare the first type to the rest of the types.
    for p ∈ rest(P, 2)
      b = r[p, :box]
      lp = local_ports[p]
      symbol_name = vars[b][lp]
      var = only(incident(decs[b], symbol_name, :name))
      type = decs[b][var, :type]
      # TODO: We use == here because we assume a type is a Symbol, like
      # :Form0, :Form1. Will a type ever be something we should compare using
      # isequal instead?
      type == type₁ || let
        var_name =  decs[b ][var, :name]
        var_name₁ = decs[b₁][var, :name]
        decapode_name =  r[b,  :name]
        decapode_name₁ = r[b₁, :name]
        error("The type of $(var_name), $(type), in decapode "*
          "\"$(decapode_name)\" does not match the type of $(var_name₁), "*
          "$(type₁), in decapode \"$(decapode_name₁)\". "*
          "(Also, check that the order of the decapodes you supplied "*
          "matches the the order you specified in the relation.)")
      end
    end
  end
end

OpenNamedDecapodeOb, OpenNamedDecapode = OpenACSetTypes(NamedDecapode, :Var)

# TODO: This does not work:
# function OpenNamedDecapode(relation, decapode, box)
function MakeOpenNamedDecapode(relation::RelationDiagram,
  decapode::NamedDecapode{Any, Any, Symbol}, box)
  P = ports(relation, box)
  J = relation[P, :junction]
  V = relation[J, :variable]
  FinFunctions = map(V) do v
    FinFunction(incident(decapode, v, :name), nparts(decapode, :Var))
  end
  OpenNamedDecapode{Any, Any, Symbol}(decapode, FinFunctions...)
end

"""    function oapply(relation::RelationDiagram, decapodes_vars::Vector{Tuple{NamedDecapode{Any, Any, Symbol}, Vector{Symbol}}})

Compose a list of decapodes as specified by the given relation diagram.

The decapodes must be given in the same order as they were specified in the
relation.

State variables (such as the (C,V) given in the head of the following
@relation) do not affect the result of a composition.

# Examples
```julia-repl
julia> compose_diff_adv = @relation (C,V) begin
  diffusion(C, ϕ₁)
  advection(C, ϕ₂, V)
  superposition(ϕ₁, ϕ₂, ϕ, C)
end;

julia> oapply(compose_diff_adv, [(Diffusion, [:C, :ϕ]),
  (Advection, [:C, :ϕ, :V]), (Superposition, [:ϕ₁, :ϕ₂, :ϕ, :C])]);
```
"""
function oapply(relation::RelationDiagram, decapodes_vars::Vector{
  Tuple{NamedDecapode{Any, Any, Symbol}, Vector{Symbol}}})
  r = relation
  copies = @. copy(first(decapodes_vars))

  # Check that the number of decapodes given matches the number of boxes in the
  # relation.
  num_boxes = nboxes(r)
  num_decapodes = length(decapodes_vars)
  # TODO: Should this be an ArgumentError?
  num_boxes == num_decapodes || error(
    "$(num_boxes) decapodes were specified in the relation but only "*
    "$(num_decapodes) were given.")

  # Check that the number of variables given in the relation is the same as the
  # number of symbols in the corresponding vector of Vars.
  # TODO: Should this be an ArgumentError?
  for b ∈ boxes(r)
    # Note: This only returns the first length mismatch found.
    num_junctions = length(incident(r, b, :box))
    num_symbols = length(decapodes_vars[b][2])
    num_junctions == num_symbols || let decapode_name = r[b,  :name]
      error("Decapode \"$(decapode_name)\" needs $(num_junctions) "*
      "variables, but $(num_symbols) were specified.")
    end
  end

  # Determine the mapping of global ports to local ports. In a RelationDiagram,
  # this is baked into the order of rows in the Port table.
  # This is a column that one could hcat to the Port table.
  local_ports = [lp for b=boxes(r) for lp=eachindex(ports(r, b))]

  # Check that types of variables connected by the same junction match.
  type_check_decapodes_composition(decapodes_vars, relation, local_ports)

  # Do namespacing.
  # Append each Var name with the name @relation gave the decapode.
  for b ∈ boxes(r)
    box_name = r[b, :name]
    for v ∈ parts(copies[b], :Var)
      var_name = copies[b][v, :name]
      copies[b][v, :name] = Symbol(box_name, '_', var_name)
    end
  end

  # Write over the name fields to be what was specified by @relation. (oapply
  # cannot combine objects whose attributes are not equal.)
  for p ∈ ports(r)
    b = r[p, :box]
    j = r[p, :junction]
    lp = local_ports[p]
    # Get the index of the row with this name in the Var.
    symbol_name = decapodes_vars[b][2][lp]
    name = Symbol((r[:name][b]), '_', symbol_name)
    # Note: only is not necessary but is a useful check the decapode is
    # well-formed. If we ever want e.g. X:Form0 and X:Form1 in a single
    # decapode, this will need refactoring.
    var = only(incident(copies[b], name, :name))
    copies[b][var, :name] = r[j, :variable]
  end

  # Compose
  apex(oapply(relation, map(boxes(r)) do b
    MakeOpenNamedDecapode(r, copies[b], b)
  end))
end

function oapply(relation::RelationDiagram,
  decapode_vars::Tuple{NamedDecapode{Any, Any, Symbol}, Vector{Symbol}})
  oapply(relation, [decapode_vars])
end

