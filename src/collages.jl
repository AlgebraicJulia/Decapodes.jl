
abstract type AbstractDecapodeMorphism end

struct BCMorphism <: AbstractDecapodeMorphism
  morphism::ACSetTransformation
end

struct ICMorphism <: AbstractDecapodeMorphism
  morphism::ACSetTransformation
end

abstract type AbstractCollage end

struct Collage <: AbstractCollage
  bc::BCMorphism
  ic::ICMorphism
end

"""    function tranfer_targets!(d::SummationDecapode, from::Int, to::Int)

Transfer any operations of which the Var `from` is the target of to `to`, excluding the special `∂ₜ` operator.
"""
function transfer_targets!(d::SummationDecapode, from::Int, to::Int)
  # All Σ of which `from` is the target:
  summations = incident(d, from, :sum)
  # All Op1s of which `from` is the target, excluding ∂ₜ:
  op1s = filter(i -> d[i, :name] != :∂ₜ, incident(d, from, :tgt))
  # All Op2s of which `from` is the target:
  op2s = incident(d, from, :res)

  # Transfer sums:
  # TODO: Change this fill idiom when ACSet assignment syntax changes.
  d[summations, :sum] = fill(to, length(summations))
  # Transfer op1s:
  d[op1s, :tgt] = fill(to, length(op1s))
  # Transfer op2s:
  d[op2s, :res] = fill(to, length(op2s))

  d
end

"""    function collate(dm::BCMorphism)

"Compile" a collage of Decapodes to a simulatable one.
```
"""
function collate(dm::BCMorphism)
  dm = dm.morphism
  d = SummationDecapode{Any, Any, Symbol}()
  copy_parts!(d, dm.codom, (:Var, :TVar, :Op1, :Op2, :Σ, :Summand))

  for (i,x) in enumerate(dm.components.Var.func)
    op_name = Symbol("∂_mask")
    mask_var = add_part!(d, :Var, type=dm.dom[i, :type], name=dm.dom[i, :name])
    tgt_name = dm.codom[x, :name]
    tgt_idx = only(incident(d, tgt_name, :name))
    d[tgt_idx, :name] = Symbol(string(d[tgt_idx, :name]) * string(i))
    res_var = add_part!(d, :Var, type=dm.codom[x, :type], name=tgt_name)
    is_tvar = !isempty(incident(d, tgt_idx, :incl))
    @show tgt_idx, mask_var, res_var
    if is_tvar
      add_part!(d, :Op2, proj1=tgt_idx, proj2=mask_var, res=res_var, op2=op_name)
    else
      transfer_targets!(d, tgt_idx, res_var)
      add_part!(d, :Op2, proj1=res_var, proj2=mask_var, res=tgt_idx, op2=op_name)
    end
    tan_dir = is_tvar ? :tgt : :src
    tangent_op1s = filter(x -> d[x, :op1]==:∂ₜ, incident(d, tgt_idx, tan_dir))
    isempty(tangent_op1s) && continue
    d[only(tangent_op1s), tan_dir] = res_var
    d[incident(d, tgt_idx, :incl), :incl] = res_var
  end
  d
end

# TODO: Make changes to the function that takes `generate` such that it defines
# the mask application function `∂_mask` by default.
"""    function make_bc_loader(dm::Collage, dimension)

Given a collage, return a function that accepts a mapping of Vars to masking functions.

This returned function will take in a mesh and non-boundary constants-and-parameters, and return a named tuple containing the parameters as evaluated on the mesh, and the regular parameters. (This final named tuple is suitable to be passed to an ODEProblem as `p`.)
"""
function make_bc_loader(dm::Collage, dimension)
  function loader(mask_funcs::Dict{Symbol, N}) where {N <: Function}
    function generator(sd, cs_ps::NamedTuple)
      vars = keys(mask_funcs)
      for (_, bc_var) in enumerate(dm.bc.morphism.dom[:name])
        bc_var ∉ vars && error("BC Variable $(string(bc_var)) is not given a generating function.")
        # TODO: Also perform error checking on the return type of the
        # corresponding mask_func. i.e. Base.return_types
        # Parameter return types should be <: Function, and Constants should be
        # tuples. Base.return_types behavior is not guaranteed to not return
        # Any, though.
      end
      mask_pairs = map(values(mask_funcs)) do func
        func(sd)
      end
      merge(
        NamedTuple{Tuple(vars)}(mask_pairs),
        cs_ps)
    end
  end
  loader
end

"""    function make_ic_loader(dm::Collage, dimension)

Given a collage, return a function that accepts a mapping of Vars to initial conditions.

This returned function will take in a mesh, and return a MultScaleArray as evaluated on the mesh. (This final MultiScaleArray is suitable to be passed to an ODEProblem as `u`.)
"""
function make_ic_loader(dm::Collage, dimension)
  function loader(ic_funcs::Dict{Symbol, N}) where {N <: Function}
    function generator(sd)
      vars = keys(ic_funcs)
      for ic_var in dm.ic.morphism.dom[:name]
        ic_var ∉ vars && error("IC Variable $(string(ic_var)) is not given a generating function.")
      end
      ics = map(values(ic_funcs)) do func
        func(sd)
      end
      codom_names = map(collect(vars)) do var
        var_idx = incident(dm.ic.morphism.dom, var, :name)
        only(dm.ic.morphism.codom[dm.ic.morphism.components.Var.func[var_idx], :name])
      end
      for (var, ic) in zip(vars,ics)
        var_idx = incident(dm.ic.morphism.dom, var, :name)
        type = only(dm.ic.morphism.dom[var_idx, :type])
        simplex = Decapodes.form_simplex(type, dimension)
        if simplex != :AllocVecCall_Error &&
          (dm.ic.morphism.dom[var_idx, :type] ∉ [:Constant, :Parameter, :infer]) &&
          length(ic) != nparts(sd, simplex)
            # TODO: This error-catch may not work if say, the number of edges is
            # equal to the number of vertices, (i.e. in a circle).
            error("IC Variable $(string(var)) was declared to be a $(string(type)), but is of length $(length(ic)), not $(nparts(sd, simplex)).")
        end
      end
      construct(PhysicsState,
        VectorForm.(ics),
        Float64[],
        collect(codom_names))
    end
  end
  loader
end

"""    function simulation_helper(dm::Collage, dimension=2)

Given a collage, return functions to load boundary conditions, initial conditions, and a simulation.
```
"""
simulation_helper(dm::Collage; dimension=2) = (
  make_bc_loader(dm,    dimension),
  make_ic_loader(dm,    dimension),
  gensim(        dm.bc, dimension=dimension))
