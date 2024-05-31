module StructuredNamedTuples

using SciMLStructures: Tunable, Constants, Caches, Discrete
import SciMLStructures: canonicalize, isscimlstructure

struct TunableWrapper
  vals
end
struct ConstantsWrapper
  vals
end
struct CachesWrapper
  vals
end
struct DiscreteWrapper
  vals
end

isscimlstructure(::NamedTuple) = true

function canonicalize(::Tunable, p::NamedTuple)
  values = map(x -> x.vals, filter(Tuple(p)) do v
    typeof(v) == TunableWrapper
  end)
  values = map(x -> x.vals, Tuple(keys_values))
  lens = length.(values)
  # TODO: Use metaprogramming instead of hardcoding AbstractVector.
  function repack(new_values::AbstractVector)
    lens_idx = 1
    buf_idx = 1
    map(keys(p), p) do k,v
      if typeof(v) == TunableWrapper
        new_value = new_values[buf_idx:buf_idx+lens[lens_idx]-1]
        buf_idx += lens[lens_idx]
        lens_idx += 1
        k,TunableWrapper(new_value)
      else
        k,v
      end
    end |> NamedTuple
  end
  return reduce(vcat, values), repack, true
end

# TODO: Use metaprogramming instead of copying-and-pasting
function canonicalize(::Constants, p::NamedTuple)
  keys_values = filter(keys(p), p) do k,v
    typeof(v) == ConstantsWrapper
  end
  values = map(x -> x.vals, Tuple(keys_values))
  lens = length.(values)
  # TODO: Use metaprogramming instead of hardcoding AbstractVector.
  function repack(new_values::AbstractVector)
    lens_idx = 1
    buf_idx = 1
    map(keys(p), p) do k,v
      if typeof(v) == ConstantsWrapper
        new_value = new_values[buf_idx:buf_idx+lens[lens_idx]-1]
        buf_idx += lens[lens_idx]
        lens_idx += 1
        k,ConstantsWrapper(new_value)
      else
        k,v
      end
    end |> NamedTuple
  end
  return reduce(vcat, values), repack, true
end

function canonicalize(::Caches, p::NamedTuple)
  keys_values = filter(keys(p), p) do k,v
    typeof(v) == CachesWrapper
  end
  values = map(x -> x.vals, Tuple(keys_values))
  lens = length.(values)
  # TODO: Use metaprogramming instead of hardcoding AbstractVector.
  function repack(new_values::AbstractVector)
    lens_idx = 1
    buf_idx = 1
    map(keys(p), p) do k,v
      if typeof(v) == CachesWrapper
        new_value = new_values[buf_idx:buf_idx+lens[lens_idx]-1]
        buf_idx += lens[lens_idx]
        lens_idx += 1
        k,CachesWrapper(new_value)
      else
        k,v
      end
    end |> NamedTuple
  end
  return reduce(vcat, values), repack, true
end

function canonicalize(::Discrete, p::NamedTuple)
  keys_values = filter(keys(p), p) do k,v
    typeof(v) == DiscreteWrapper
  end
  values = map(x -> x.vals, Tuple(keys_values))
  lens = length.(values)
  # TODO: Use metaprogramming instead of hardcoding AbstractVector.
  function repack(new_values::AbstractVector)
    lens_idx = 1
    buf_idx = 1
    map(keys(p), p) do k,v
      if typeof(v) == DiscreteWrapper
        new_value = new_values[buf_idx:buf_idx+lens[lens_idx]-1]
        buf_idx += lens[lens_idx]
        lens_idx += 1
        k,DiscreteWrapper(new_value)
      else
        k,v
      end
    end |> NamedTuple
  end
  return reduce(vcat, values), repack, true
end

end
