
import Catlab.WiringDiagrams: oapply

OpenNamedDecapodeOb, OpenNamedDecapode = OpenACSetTypes(NamedDecapode, :Var)

#FIXME: why can't we just add a constructor for OpenNamedDecapode
function OpenPode(d::NamedDecapode, names::AbstractVector{Symbol})
  legs = map(names) do name
    FinFunction(incident(d, name, :name), nparts(d, :Var))
  end
  OpenNamedDecapode{Any, Any, Symbol}(d, legs...)
end

apex(decapode::OpenNamedDecapode) = apex(decapode.cospan)
legs(decapode::OpenNamedDecapode) = legs(decapode.cospan)
feet(decapode::OpenNamedDecapode) = decapode.feet

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


# TODO: This does not work:
# function OpenNamedDecapode(relation, decapode, box)
function MakeOpenNamedDecapode(relation::RelationDiagram, decapode::NamedDecapode, box)
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
function oapply(relation::RelationDiagram, decapodes_vars::Vector{OpenNamedDecapode})
  r = relation
  # FIXME: in this line, you should cast the NamedDecapode{S,T, Symbol} to NamedDecapode{S,T,Vector{Symbol}}
  # This will allow you to return namespace scoped variables.
  copies = deepcopy(decapodes_vars)

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
    num_symbols = length(feet(decapodes_vars[b]))
    num_junctions == num_symbols || let decapode_name = r[b,  :name]
      error("Component $(decapode_name) expects $(num_junctions) interface variables, but number of feet is $(num_symbols).")
    end
    # FIXME: this should also check the types of the feet.
  end

  # Determine the mapping of global ports to local ports. In a RelationDiagram,
  # this is baked into the order of rows in the Port table.
  # This is a column that one could hcat to the Port table.
  local_ports = [lp for b=boxes(r) for lp=eachindex(ports(r, b))]

  # Check that types of variables connected by the same junction match.
  type_check_decapodes_composition(decapodes_vars, relation, local_ports)

  # Do namespacing.
  # Append each Var name with the name @relation gave the decapode.
  # FIXME: return a list of symbols here instead of underscore separated list.
  for b ∈ boxes(r)
    box_name = r[b, :name]
    for v ∈ parts(copies[b], :Var)
      var_name = copies[b][v, :name]
      copies[b][v, :name] = Symbol(box_name, '_', var_name)
    end
  end

  # Write over the name fields to be what was specified by @relation. (oapply
  # cannot combine objects whose attributes are not equal.)

  newnames = map(boxes(r)) do b
    pode = decapodes_vars[b]
    ports = incident(r, b, :box)
    localnames = map(enumerate(ports)) do (i,p)
      localnamevec = feet(pode)[i][:name]
      length(localnamevec) == 1 || error("Decapode arguments to oapply do not support bundling. Each foot should have 1 vertex")
      localnamevec[1]
    end

    map(zip(ports, localnames)) do (p, lname)
      # FIXME: this will break when we add proper namespacing
      # Note: only is not necessary but is a useful check the decapode is
      # well-formed. If we ever want e.g. X:Form0 and X:Form1 in a single
      # decapode, this will need refactoring.
      name = Symbol(r[b, :name], '_', lname)
      var = only(incident(pode, name, :name))
      j = r[p, :junction]
      globalname = r[j, :variable]
      pode[var, :name] = globalname
      return globalname
    end
  end

  newpodes = map(boxes(r)) do b
    OpenPode(apex(decapodes_vars[b]), newnames[b])
  end

  # Compose
  oapply(relation, newpodes)
end

# function oapply(relation::RelationDiagram,
#   decapode_vars::Tuple{NamedDecapode{Any, Any, Symbol}, Vector{Symbol}})
#   oapply(relation, [decapode_vars])
# end

