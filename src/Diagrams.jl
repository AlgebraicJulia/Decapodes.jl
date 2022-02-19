module Diagrams

using CombinatorialSpaces.ExteriorCalculus
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.Programs.DiagrammaticPrograms: parse_diagram, parse_gat_expr
using Unicode

export @decapode, Decapodes2D

macro decapode(cat, body)
  :(parse_exp_diagram($(esc(cat)), $(Meta.quot(body)), free=true))
end

function parse_exp_diagram(cat, body; kw...)
  objs = filter(exp -> !(exp isa LineNumberNode) && (exp.head == :(::)), body.args)
  obj_map = Dict{Symbol, Expr}()
  for obj in objs
    names = obj.args[1]
    type = obj.args[2]
    # Expand if object is a tuple or single object
    if names isa Symbol
      obj_map[names] = type
    else
      for name in names.args
        obj_map[name] = type
      end
    end
  end
  new_body = quote end
  for exp in body.args
    if !(exp isa LineNumberNode) && exp.head == :call && exp.args[1] == :(==)
      exprs = []
      if exp.args[2] isa Symbol
        res, _ = expand_expr!(exprs, exp.args[3], cat, obj_map; make_var = false)
        push!(exprs, :($(exp.args[2]) == $(res)))
      elseif exp.args[3] isa Symbol
        res, _ = expand_expr!(exprs, exp.args[2], cat, obj_map; make_var = false)
        push!(exprs, :($(exp.args[3]) == $(res)))
      else
        res1, _ = expand_expr!(exprs, exp.args[2], cat, obj_map; make_var = true)
        res2, _ = expand_expr!(exprs, exp.args[3], cat, obj_map; make_var = false)
        push!(exprs, :($(res1) == $(res2)))
      end
      append!(new_body.args, exprs)
    else
      push!(new_body.args, exp)
    end
  end
  parse_diagram(cat, new_body; kw...)
end

i2sub = Dict(
  '0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', '5'=>'₅',
  '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉', '-'=>'₋'
)
i2sup = Dict(
  '0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', '5'=>'⁵',
  '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻'
)

form_sup(type::Expr) = form_sup(type.args[1])
form_sub(type::Expr) = form_sub(type.args[1])

form_sup(type::Symbol) = startswith("$type", "Dual") ?
                                               i2sup[last("$type")] * '̃' :  i2sup[last("$type")]
form_sub(type::Symbol) = startswith("$type", "Dual") ? i2sub[last("$type")] * '̃' : i2sub[last("$type")]


# Currently every projection needs its own index
function make_proj(ind, types)
  op = "proj" * i2sub["$ind"[1]] * "_" * join(form_sup.(types)) * form_sub(types[ind])
  Symbol(Unicode.normalize(op))
end

function parse_type(type)
  if head(type) == :otimes
    Meta.parse("otimes{$(join(["$(head(t)){$(t.args[1])}" for t in type.args], ","))}")
  else
    Meta.parse("$(head(type)){$(type.args[1])}")
  end
end


get_name(expr) = expr isa Symbol ? expr : get_name(expr.args[1])

function expand_expr!(expr_arr, expr, cat, obj_map; make_var = true)
  if expr isa Symbol
    expr, obj_map[expr]
  elseif length(expr.args) > 2
    if expr.args[1] == :+
      # Only process 2 of these at a time (sum first, evaluate second as sum)
      res1 = expand_expr!(expr_arr, expr.args[2], cat, obj_map; make_var = true)

      res2 = length(expr.args) > 3 ?
      expand_expr!(expr_arr, :(+($(expr.args[3:end]...))), cat, obj_map; make_var = true) :
              expand_expr!(expr_arr, expr.args[3], cat, obj_map; make_var = true)
      type = res1[2]

      sum_name = Symbol(Unicode.normalize("sum" * form_sub(type)))
      new_var = gensym()
      nv_type = :(otimes{$(res1[2]), $(res2[2])})
      push!(expr_arr, :($(new_var)::$(nv_type)))
      push!(expr_arr, :($(res1[1]) == $(make_proj(1, [res1[2], res2[2]]))($(new_var))))
      push!(expr_arr, :($(res2[1]) == $(make_proj(2, [res1[2], res2[2]]))($(new_var))))
      if make_var
        name = gensym()
        push!(expr_arr, :($(name)::$(type)))
        push!(expr_arr, :($(name) == $(sum_name)($(new_var))))
        name, type
      else
        :($(sum_name)($(new_var))), type
      end
    else
      res = map(expr.args[2:end]) do ex
        expand_expr!(expr_arr, ex, cat, obj_map; make_var = true)
      end
      new_var = gensym()
      nv_type = Meta.parse("otimes{$(join(last.(res), ","))}")
      push!(expr_arr, :($(new_var)::$(nv_type)))
      append!(expr_arr, map(enumerate(res)) do (i, obj)
                :($(obj[1]) == $(make_proj(i, last.(res)))($(new_var)))
              end)
      type = parse_type(codom(parse_gat_expr(FinCat(cat), expr.args[1])))
      if make_var
        name = gensym()
        push!(expr_arr, :($(name)::$(type)))
        push!(expr_arr, :($(name) == $(expr.args[1])($(new_var))))
        name, type
      else
        :($(expr.args[1])($(new_var))), type
      end
    end
  else
    type = parse_type(codom(parse_gat_expr(FinCat(cat), expr.args[1])))
    arg_var, arg_type = expand_expr!(expr_arr, expr.args[2], cat, obj_map; make_var = false)
    expr.args[2] = arg_var
    cur_var = expr
    if make_var
      name = gensym()
      push!(expr_arr, :($(name)::$(type)))
      push!(expr_arr, :($(name) == $(expr)))
      cur_var = name
    end
    cur_var, type
  end
end

""" Decapodes2D
A schema which includes any homomorphisms that may be added by the @decapode
macro.

TODO: This should be chipped away at as more of this tooling takes advantage
of the Catlab GAT system
"""
@present Decapodes2D(FreeExtCalc2D) begin
  X::Space
  proj₁_⁰⁰₀::Hom(Form0(X)⊗Form0(X),Form0(X))
  proj₂_⁰⁰₀::Hom(Form0(X)⊗Form0(X),Form0(X))
  proj₁_⁰⁰₀⁺::Hom(Form0(X)⊕Form0(X),Form0(X))
  proj₂_⁰⁰₀⁺::Hom(Form0(X)⊕Form0(X),Form0(X))
  proj₁_⁰¹₀::Hom(Form0(X)⊗Form1(X),Form0(X))
  proj₂_⁰¹₁::Hom(Form0(X)⊗Form1(X),Form1(X))
  proj₁_⁰¹₀⁺::Hom(Form0(X)⊕Form1(X),Form0(X))
  proj₂_⁰¹₁⁺::Hom(Form0(X)⊕Form1(X),Form1(X))
  proj₁_⁰²₀::Hom(Form0(X)⊗Form2(X),Form0(X))
  proj₂_⁰²₂::Hom(Form0(X)⊗Form2(X),Form2(X))
  proj₁_⁰²₀⁺::Hom(Form0(X)⊕Form2(X),Form0(X))
  proj₂_⁰²₂⁺::Hom(Form0(X)⊕Form2(X),Form2(X))
  proj₁_⁰⁰̃₀::Hom(Form0(X)⊗DualForm0(X),Form0(X))
  proj₂_⁰⁰̃₀̃::Hom(Form0(X)⊗DualForm0(X),DualForm0(X))
  proj₁_⁰⁰̃₀⁺::Hom(Form0(X)⊕DualForm0(X),Form0(X))
  proj₂_⁰⁰̃₀̃⁺::Hom(Form0(X)⊕DualForm0(X),DualForm0(X))
  proj₁_⁰¹̃₀::Hom(Form0(X)⊗DualForm1(X),Form0(X))
  proj₂_⁰¹̃₁̃::Hom(Form0(X)⊗DualForm1(X),DualForm1(X))
  proj₁_⁰¹̃₀⁺::Hom(Form0(X)⊕DualForm1(X),Form0(X))
  proj₂_⁰¹̃₁̃⁺::Hom(Form0(X)⊕DualForm1(X),DualForm1(X))
  proj₁_⁰²̃₀::Hom(Form0(X)⊗DualForm2(X),Form0(X))
  proj₂_⁰²̃₂̃::Hom(Form0(X)⊗DualForm2(X),DualForm2(X))
  proj₁_⁰²̃₀⁺::Hom(Form0(X)⊕DualForm2(X),Form0(X))
  proj₂_⁰²̃₂̃⁺::Hom(Form0(X)⊕DualForm2(X),DualForm2(X))
  proj₁_¹⁰₁::Hom(Form1(X)⊗Form0(X),Form1(X))
  proj₂_¹⁰₀::Hom(Form1(X)⊗Form0(X),Form0(X))
  proj₁_¹⁰₁⁺::Hom(Form1(X)⊕Form0(X),Form1(X))
  proj₂_¹⁰₀⁺::Hom(Form1(X)⊕Form0(X),Form0(X))
  proj₁_¹¹₁::Hom(Form1(X)⊗Form1(X),Form1(X))
  proj₂_¹¹₁::Hom(Form1(X)⊗Form1(X),Form1(X))
  proj₁_¹¹₁⁺::Hom(Form1(X)⊕Form1(X),Form1(X))
  proj₂_¹¹₁⁺::Hom(Form1(X)⊕Form1(X),Form1(X))
  proj₁_¹²₁::Hom(Form1(X)⊗Form2(X),Form1(X))
  proj₂_¹²₂::Hom(Form1(X)⊗Form2(X),Form2(X))
  proj₁_¹²₁⁺::Hom(Form1(X)⊕Form2(X),Form1(X))
  proj₂_¹²₂⁺::Hom(Form1(X)⊕Form2(X),Form2(X))
  proj₁_¹⁰̃₁::Hom(Form1(X)⊗DualForm0(X),Form1(X))
  proj₂_¹⁰̃₀̃::Hom(Form1(X)⊗DualForm0(X),DualForm0(X))
  proj₁_¹⁰̃₁⁺::Hom(Form1(X)⊕DualForm0(X),Form1(X))
  proj₂_¹⁰̃₀̃⁺::Hom(Form1(X)⊕DualForm0(X),DualForm0(X))
  proj₁_¹¹̃₁::Hom(Form1(X)⊗DualForm1(X),Form1(X))
  proj₂_¹¹̃₁̃::Hom(Form1(X)⊗DualForm1(X),DualForm1(X))
  proj₁_¹¹̃₁⁺::Hom(Form1(X)⊕DualForm1(X),Form1(X))
  proj₂_¹¹̃₁̃⁺::Hom(Form1(X)⊕DualForm1(X),DualForm1(X))
  proj₁_¹²̃₁::Hom(Form1(X)⊗DualForm2(X),Form1(X))
  proj₂_¹²̃₂̃::Hom(Form1(X)⊗DualForm2(X),DualForm2(X))
  proj₁_¹²̃₁⁺::Hom(Form1(X)⊕DualForm2(X),Form1(X))
  proj₂_¹²̃₂̃⁺::Hom(Form1(X)⊕DualForm2(X),DualForm2(X))
  proj₁_²⁰₂::Hom(Form2(X)⊗Form0(X),Form2(X))
  proj₂_²⁰₀::Hom(Form2(X)⊗Form0(X),Form0(X))
  proj₁_²⁰₂⁺::Hom(Form2(X)⊕Form0(X),Form2(X))
  proj₂_²⁰₀⁺::Hom(Form2(X)⊕Form0(X),Form0(X))
  proj₁_²¹₂::Hom(Form2(X)⊗Form1(X),Form2(X))
  proj₂_²¹₁::Hom(Form2(X)⊗Form1(X),Form1(X))
  proj₁_²¹₂⁺::Hom(Form2(X)⊕Form1(X),Form2(X))
  proj₂_²¹₁⁺::Hom(Form2(X)⊕Form1(X),Form1(X))
  proj₁_²²₂::Hom(Form2(X)⊗Form2(X),Form2(X))
  proj₂_²²₂::Hom(Form2(X)⊗Form2(X),Form2(X))
  proj₁_²²₂⁺::Hom(Form2(X)⊕Form2(X),Form2(X))
  proj₂_²²₂⁺::Hom(Form2(X)⊕Form2(X),Form2(X))
  proj₁_²⁰̃₂::Hom(Form2(X)⊗DualForm0(X),Form2(X))
  proj₂_²⁰̃₀̃::Hom(Form2(X)⊗DualForm0(X),DualForm0(X))
  proj₁_²⁰̃₂⁺::Hom(Form2(X)⊕DualForm0(X),Form2(X))
  proj₂_²⁰̃₀̃⁺::Hom(Form2(X)⊕DualForm0(X),DualForm0(X))
  proj₁_²¹̃₂::Hom(Form2(X)⊗DualForm1(X),Form2(X))
  proj₂_²¹̃₁̃::Hom(Form2(X)⊗DualForm1(X),DualForm1(X))
  proj₁_²¹̃₂⁺::Hom(Form2(X)⊕DualForm1(X),Form2(X))
  proj₂_²¹̃₁̃⁺::Hom(Form2(X)⊕DualForm1(X),DualForm1(X))
  proj₁_²²̃₂::Hom(Form2(X)⊗DualForm2(X),Form2(X))
  proj₂_²²̃₂̃::Hom(Form2(X)⊗DualForm2(X),DualForm2(X))
  proj₁_²²̃₂⁺::Hom(Form2(X)⊕DualForm2(X),Form2(X))
  proj₂_²²̃₂̃⁺::Hom(Form2(X)⊕DualForm2(X),DualForm2(X))
  proj₁_⁰̃⁰₀̃::Hom(DualForm0(X)⊗Form0(X),DualForm0(X))
  proj₂_⁰̃⁰₀::Hom(DualForm0(X)⊗Form0(X),Form0(X))
  proj₁_⁰̃⁰₀̃⁺::Hom(DualForm0(X)⊕Form0(X),DualForm0(X))
  proj₂_⁰̃⁰₀⁺::Hom(DualForm0(X)⊕Form0(X),Form0(X))
  proj₁_⁰̃¹₀̃::Hom(DualForm0(X)⊗Form1(X),DualForm0(X))
  proj₂_⁰̃¹₁::Hom(DualForm0(X)⊗Form1(X),Form1(X))
  proj₁_⁰̃¹₀̃⁺::Hom(DualForm0(X)⊕Form1(X),DualForm0(X))
  proj₂_⁰̃¹₁⁺::Hom(DualForm0(X)⊕Form1(X),Form1(X))
  proj₁_⁰̃²₀̃::Hom(DualForm0(X)⊗Form2(X),DualForm0(X))
  proj₂_⁰̃²₂::Hom(DualForm0(X)⊗Form2(X),Form2(X))
  proj₁_⁰̃²₀̃⁺::Hom(DualForm0(X)⊕Form2(X),DualForm0(X))
  proj₂_⁰̃²₂⁺::Hom(DualForm0(X)⊕Form2(X),Form2(X))
  proj₁_⁰̃⁰̃₀̃::Hom(DualForm0(X)⊗DualForm0(X),DualForm0(X))
  proj₂_⁰̃⁰̃₀̃::Hom(DualForm0(X)⊗DualForm0(X),DualForm0(X))
  proj₁_⁰̃⁰̃₀̃⁺::Hom(DualForm0(X)⊕DualForm0(X),DualForm0(X))
  proj₂_⁰̃⁰̃₀̃⁺::Hom(DualForm0(X)⊕DualForm0(X),DualForm0(X))
  proj₁_⁰̃¹̃₀̃::Hom(DualForm0(X)⊗DualForm1(X),DualForm0(X))
  proj₂_⁰̃¹̃₁̃::Hom(DualForm0(X)⊗DualForm1(X),DualForm1(X))
  proj₁_⁰̃¹̃₀̃⁺::Hom(DualForm0(X)⊕DualForm1(X),DualForm0(X))
  proj₂_⁰̃¹̃₁̃⁺::Hom(DualForm0(X)⊕DualForm1(X),DualForm1(X))
  proj₁_⁰̃²̃₀̃::Hom(DualForm0(X)⊗DualForm2(X),DualForm0(X))
  proj₂_⁰̃²̃₂̃::Hom(DualForm0(X)⊗DualForm2(X),DualForm2(X))
  proj₁_⁰̃²̃₀̃⁺::Hom(DualForm0(X)⊕DualForm2(X),DualForm0(X))
  proj₂_⁰̃²̃₂̃⁺::Hom(DualForm0(X)⊕DualForm2(X),DualForm2(X))
  proj₁_¹̃⁰₁̃::Hom(DualForm1(X)⊗Form0(X),DualForm1(X))
  proj₂_¹̃⁰₀::Hom(DualForm1(X)⊗Form0(X),Form0(X))
  proj₁_¹̃⁰₁̃⁺::Hom(DualForm1(X)⊕Form0(X),DualForm1(X))
  proj₂_¹̃⁰₀⁺::Hom(DualForm1(X)⊕Form0(X),Form0(X))
  proj₁_¹̃¹₁̃::Hom(DualForm1(X)⊗Form1(X),DualForm1(X))
  proj₂_¹̃¹₁::Hom(DualForm1(X)⊗Form1(X),Form1(X))
  proj₁_¹̃¹₁̃⁺::Hom(DualForm1(X)⊕Form1(X),DualForm1(X))
  proj₂_¹̃¹₁⁺::Hom(DualForm1(X)⊕Form1(X),Form1(X))
  proj₁_¹̃²₁̃::Hom(DualForm1(X)⊗Form2(X),DualForm1(X))
  proj₂_¹̃²₂::Hom(DualForm1(X)⊗Form2(X),Form2(X))
  proj₁_¹̃²₁̃⁺::Hom(DualForm1(X)⊕Form2(X),DualForm1(X))
  proj₂_¹̃²₂⁺::Hom(DualForm1(X)⊕Form2(X),Form2(X))
  proj₁_¹̃⁰̃₁̃::Hom(DualForm1(X)⊗DualForm0(X),DualForm1(X))
  proj₂_¹̃⁰̃₀̃::Hom(DualForm1(X)⊗DualForm0(X),DualForm0(X))
  proj₁_¹̃⁰̃₁̃⁺::Hom(DualForm1(X)⊕DualForm0(X),DualForm1(X))
  proj₂_¹̃⁰̃₀̃⁺::Hom(DualForm1(X)⊕DualForm0(X),DualForm0(X))
  proj₁_¹̃¹̃₁̃::Hom(DualForm1(X)⊗DualForm1(X),DualForm1(X))
  proj₂_¹̃¹̃₁̃::Hom(DualForm1(X)⊗DualForm1(X),DualForm1(X))
  proj₁_¹̃¹̃₁̃⁺::Hom(DualForm1(X)⊕DualForm1(X),DualForm1(X))
  proj₂_¹̃¹̃₁̃⁺::Hom(DualForm1(X)⊕DualForm1(X),DualForm1(X))
  proj₁_¹̃²̃₁̃::Hom(DualForm1(X)⊗DualForm2(X),DualForm1(X))
  proj₂_¹̃²̃₂̃::Hom(DualForm1(X)⊗DualForm2(X),DualForm2(X))
  proj₁_¹̃²̃₁̃⁺::Hom(DualForm1(X)⊕DualForm2(X),DualForm1(X))
  proj₂_¹̃²̃₂̃⁺::Hom(DualForm1(X)⊕DualForm2(X),DualForm2(X))
  proj₁_²̃⁰₂̃::Hom(DualForm2(X)⊗Form0(X),DualForm2(X))
  proj₂_²̃⁰₀::Hom(DualForm2(X)⊗Form0(X),Form0(X))
  proj₁_²̃⁰₂̃⁺::Hom(DualForm2(X)⊕Form0(X),DualForm2(X))
  proj₂_²̃⁰₀⁺::Hom(DualForm2(X)⊕Form0(X),Form0(X))
  proj₁_²̃¹₂̃::Hom(DualForm2(X)⊗Form1(X),DualForm2(X))
  proj₂_²̃¹₁::Hom(DualForm2(X)⊗Form1(X),Form1(X))
  proj₁_²̃¹₂̃⁺::Hom(DualForm2(X)⊕Form1(X),DualForm2(X))
  proj₂_²̃¹₁⁺::Hom(DualForm2(X)⊕Form1(X),Form1(X))
  proj₁_²̃²₂̃::Hom(DualForm2(X)⊗Form2(X),DualForm2(X))
  proj₂_²̃²₂::Hom(DualForm2(X)⊗Form2(X),Form2(X))
  proj₁_²̃²₂̃⁺::Hom(DualForm2(X)⊕Form2(X),DualForm2(X))
  proj₂_²̃²₂⁺::Hom(DualForm2(X)⊕Form2(X),Form2(X))
  proj₁_²̃⁰̃₂̃::Hom(DualForm2(X)⊗DualForm0(X),DualForm2(X))
  proj₂_²̃⁰̃₀̃::Hom(DualForm2(X)⊗DualForm0(X),DualForm0(X))
  proj₁_²̃⁰̃₂̃⁺::Hom(DualForm2(X)⊕DualForm0(X),DualForm2(X))
  proj₂_²̃⁰̃₀̃⁺::Hom(DualForm2(X)⊕DualForm0(X),DualForm0(X))
  proj₁_²̃¹̃₂̃::Hom(DualForm2(X)⊗DualForm1(X),DualForm2(X))
  proj₂_²̃¹̃₁̃::Hom(DualForm2(X)⊗DualForm1(X),DualForm1(X))
  proj₁_²̃¹̃₂̃⁺::Hom(DualForm2(X)⊕DualForm1(X),DualForm2(X))
  proj₂_²̃¹̃₁̃⁺::Hom(DualForm2(X)⊕DualForm1(X),DualForm1(X))
  proj₁_²̃²̃₂̃::Hom(DualForm2(X)⊗DualForm2(X),DualForm2(X))
  proj₂_²̃²̃₂̃::Hom(DualForm2(X)⊗DualForm2(X),DualForm2(X))
  proj₁_²̃²̃₂̃⁺::Hom(DualForm2(X)⊕DualForm2(X),DualForm2(X))
  proj₂_²̃²̃₂̃⁺::Hom(DualForm2(X)⊕DualForm2(X),DualForm2(X))
  sum₀::Hom(Form0(X)⊗Form0(X),Form0(X))
  sum₁::Hom(Form1(X)⊗Form1(X),Form1(X))
  sum₂::Hom(Form2(X)⊗Form2(X),Form2(X))
  sum₀̃::Hom(DualForm0(X)⊗DualForm0(X),DualForm0(X))
  sum₁̃::Hom(DualForm1(X)⊗DualForm1(X),DualForm1(X))
  sum₂̃::Hom(DualForm2(X)⊗DualForm2(X),DualForm2(X))
end

end
