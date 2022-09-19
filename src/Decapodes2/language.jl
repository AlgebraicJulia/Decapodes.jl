# Definition of the Decapodes DSL AST
# TODO: do functions and tans need to be parameterized by a space?
# TODO: make recursive
@data Term begin
  Var(Symbol)
  Judge(Var, Symbol, Symbol) # Symbol 1: Form0 Symbol 2: X
  Eq(Term, Term)
  AppCirc1(Vector{Symbol}, Var)
  AppCirc2(Vector{Symbol}, Var, Var)
  App1(Symbol, Var)
  App2(Symbol, Var, Var)
  Plus(Var, Var)
  Tan(Var)
end

# A struct to store a complete Decapode
# TODO: Have the decopode macro compile to a DecaExpr
struct DecaExpr
  judgements::Vector{Judge}
  equations::Vector{Eq}
end

term(s::Symbol) = Var(normalize_unicode(s))

term(expr::Expr) = begin
    @match expr begin
        Expr(a) => Var(normalize_unicode(a))
        Expr(:call, :∂ₜ, b) => Tan(term(b))
        Expr(:call, Expr(:call, :∘, a...), b) => AppCirc1(a, Var(b))
        Expr(:call, a, b) => App1(a, term(b))
        Expr(:call, Expr(:call, :∘, f...), x, y) => AppCirc1(f, Var(x), Var(y))
        Expr(:call, f, x, y) => App2(f, term(x), term(y))
        Expr(:call, :∘, a...) => (:AppCirc1, map(term, a))
        x => error("Cannot construct term from  $x")
    end
end

function parse_decapode(expr::Expr)
    stmts = map(expr.args) do line 
        @match line begin
            ::LineNumberNode => missing
            Expr(:(::), a, b) => Judge(Var(a),b.args[1], b.args[2])
            Expr(:call, :(==), lhs, rhs) => Eq(term(lhs), term(rhs))
            x => x
        end
    end |> skipmissing |> collect
    judges = []
    eqns = []
    for s in stmts
        if typeof(s) == Judge
            push!(judges, s)
        elseif typeof(s) == Eq
            push!(eqns, s)
        end
    end
    DecaExpr(judges, eqns)
end


# to_decapode helper functions
reduce_lhs!(eq::Eq, d::AbstractDecapode, syms::Dict{Symbol, Int}) =
  let ! = reduce_lhs! # This will be needed once we upgrade to a recursive grammar
    @match eq._1 begin
      Var(x) => syms[x]
      App1(f, x) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op1, src=syms[x._1], tgt=res_var, op1=f)
        return res_var
      end
      App2(f, x, y) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=res_var, op2=f)
        return res_var
      end
      AppCirc1(fs, x) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op1, src=syms[x._1], tgt=res_var, op1=fs)
        return res_var
      end
      AppCirc2(fs, x, y) => begin
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=res_var, op2=fs)
        return res_var
      end
      Tan(x) => begin
        # TODO: this is creating a spurious variablbe with the same name
        txv = add_part!(d, :Var, type=:infer)
        tx = add_part!(d, :TVar, incl=txv)
        tanop = add_part!(d, :Op1, src=syms[x._1], tgt=txv, op1=DerivOp)
        return txv #syms[x._1]
      end
      Plus(x, y) => begin # TODO: plus is an Op2 so just fold it into App2
        res_var = add_part!(d, :Var, type=:infer)
        add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=res_var, op2=:plus)
        return res_var
      end
      _ => -1 # TODO: make this throw an error or something
    end
  end
# TODO: lots of code duplication between reduce_lhs! and reduce_rhs!
# The duplicate code should be abstracted into another helper function
reduce_rhs!(eq::Eq, d::AbstractDecapode, syms::Dict{Symbol, Int}, lhs_ref::Int) =
  let ! = reduce_rhs! # Again only necessary once we upgrade language
    @match eq._2 begin
      App1(f, x) => add_part!(d, :Op1, src=syms[x._1], tgt=lhs_ref, op1=f)
      App2(f, x, y) => add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=lhs_ref, op2=f)
      AppCirc1(fs, x) => add_part!(d, :Op1, src=syms[x._1], tgt=lhs_ref, op1=fs)
      AppCirc2(fs, x, y) => add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=lhs_ref, op2=fs)
      Plus(x, y) => add_part!(d, :Op2, proj1=syms[x._1], proj2=syms[y._1], res=lhs_ref, op2=:plus)
      _ => -1 # TODO: Throw an error or handle case where RHS is a raw variable or tangent
    end
  end

""" Takes a DecaExpr (i.e. what should be constructed using the @decapode macro)
and gives a Decapode ACSet which represents equalities as two operations with the
same tgt or res map.
"""
# Just to get up and running, I tagged untyped variables with :infer
# TODO: write a type checking/inference step for this function to 
# overwrite the :infer tags
function Decapode(e::DecaExpr)
  d = Decapode{Any, Any}()
  symbol_table = Dict{Symbol, Int}()
  for judgement in e.judgements
    var_id = add_part!(d, :Var, type=(judgement._2, judgement._3))
    symbol_table[judgement._1._1] = var_id
  end
  for eq in e.equations
    v = reduce_lhs!(eq, d, symbol_table)
    reduce_rhs!(eq, d, symbol_table, v)
  end
  return d
end

function NamedDecapode(e::DecaExpr)
    d = NamedDecapode{Any, Any, Symbol}()
    symbol_table = Dict{Symbol, Int}()
    for judgement in e.judgements
      var_id = add_part!(d, :Var, name=judgement._1._1, type=judgement._2)
      symbol_table[judgement._1._1] = var_id
    end
    for eq in e.equations
      v = reduce_lhs!(eq, d, symbol_table)
      reduce_rhs!(eq, d, symbol_table, v)
    end
    fill_names!(d)
    d[:name] = map(normalize_unicode,d[:name])
    return d
end
