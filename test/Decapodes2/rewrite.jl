using Test

using Decapodes
using Catlab.Graphics
using Catlab, Catlab.Graphs, Catlab.Graphics, Catlab.CategoricalAlgebra
using Catlab.Theories, Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite
const hom = AlgebraicRewriting.homomorphism
const Var = AlgebraicRewriting.Var

draw(g; kw...) = to_graphviz(g; node_labels=true, edge_labels=true, kw...)
draw(f::ACSetTransformation; kw...) =
  to_graphviz(f; node_labels=true, edge_labels=true, draw_codom=false, kw...)


#########################

function get_target_indices(deca_source)
  targetVars = []
  for var in parts(deca_source, :Var)
    op1Count = length(incident(deca_source, var, :tgt))
    op2Count = length(incident(deca_source, var, :res))
    sumCount = length(incident(deca_source, var, :sum))

    tot = op1Count + op2Count + sumCount
    if(tot >= 2)
      append!(targetVars, var)
    end
  end

  return targetVars
end

function get_preprocess_indices(deca_source)
  targetOp2 = []
  targetSum = []

  targetVars = get_target_indices(deca_source)

  for var in targetVars
    append!(targetOp2, incident(deca_source, var, :res))
    append!(targetSum, incident(deca_source, var, :sum))
  end

  return targetOp2, targetSum
end

function preprocess_rewrite(deca_source)
  targetOp2, targetSum = get_preprocess_indices(deca_source)

  # If we don't need to preprocess then don't
  if(length(targetOp2) == 0 && length(targetSum) == 0)
    return deca_source
  end

  LHS = []
  RHS = []

  SuperMatch = []
  SuperVarMap = Vector{Int64}()
  SuperOp2Map = Vector{Int64}()

  # Process all of the target rewrites for op2
  for opID in targetOp2

    vars = [deca_source[opID, :proj1], deca_source[opID, :proj2], deca_source[opID, :res]]
    types = deca_source[vars, :type]
    names = deca_source[vars, :name]
    op2name = deca_source[opID, :op2]

    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 3
      type = types
      name = names

      Op2 = 1
      proj1 = [1]
      proj2 = [2]
      res = [3]
      op2 = op2name
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 3
      type = types
      name = names
    end

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 4
      type = vcat(types, types[end])
      name = vcat(names, :temp)

      Op1 = 1
      src = [4]
      tgt = [3]
      op1 = [:temp]

      Op2 = 1
      proj1 = [1]
      proj2 = [2]
      res = [4]
      op2 = op2name
    end

    L = ACSetTransformation(I, Match, Var = 1:3)
    R = ACSetTransformation(I, Sub, Var = 1:3)

    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    push!(SuperOp2Map, opID)
  end
    

  # Combine all rules in parallel and apply
  rule = Rule(oplus(LHS), oplus(RHS))
  m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op2 = SuperOp2Map)

  rewrite_match(rule, m)
end

function rewrite_decapode(deca_source)
  # Just for now, I'll get it working with op1 only
  # Considering op2 and summations will make this significantly more difficult

  targetVars = get_target_indices(deca_source)

  if(length(targetVars) == 0)
    return deca_source
  end

  LHS = []
  RHS = []

  SuperMatch = []
  SuperVarMap = Vector{Int64}()
  SuperOp1Map = Vector{Int64}()

  for varID in targetVars
    targetOp1 = incident(deca_source, varID, :tgt)
    vars = vcat(deca_source[targetOp1, :src], varID)

    num_nodes_match = length(vars)
    nary_of_rewrite = num_nodes_match - 1

    result_index = num_nodes_match
    sum_index = 2 * result_index

    variable_types = deca_source[vars, :type]
    variable_var =  deca_source[vars, :name]
    variable_op1 = deca_source[targetOp1, :op1]

    # Cannot work with undefined types
    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = num_nodes_match
      type = variable_types
      name = variable_var

      # This will probably break for rewrites including 
      # Non-Op1 rewrites
      Op1 = nary_of_rewrite
      src = 1:nary_of_rewrite
      tgt = fill(result_index, nary_of_rewrite)
      op1 = variable_op1
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = num_nodes_match
      type = variable_types
      name = variable_var
    end 

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 2 * num_nodes_match
      type = vcat(variable_types, variable_types)
      name = vcat(variable_var, map(x -> Symbol("••",x), 1:nary_of_rewrite), [:sum])
      Op1 = nary_of_rewrite + 1
      src = vcat(1:nary_of_rewrite, sum_index)
      tgt = vcat(num_nodes_match+1:sum_index-1, [result_index])
      op1 = vcat(variable_op1, Symbol(:avg, nary_of_rewrite))
      Σ = 1
      sum = [sum_index]
      Summand = nary_of_rewrite
      summand = num_nodes_match+1:sum_index-1
      summation = fill(1, nary_of_rewrite)
    end

    L = ACSetTransformation(I, Match, Var = 1:num_nodes_match);
    R = ACSetTransformation(I, Sub, Var = 1:num_nodes_match);
    
    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    append!(SuperOp1Map, targetOp1)
  end

  rule = Rule(oplus(LHS), oplus(RHS))

  m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op1 = SuperOp1Map)
  rewrite_match(rule, m)
end

"""
Test for multiple op1 from the same source
Test for multiple potential rewrites
Test for no valid rewrite
Test for large nary
"""
# Trivial test
# Make sure var, op names and forms are preserved
DecaTest1 = quote
  D₁::Form1{X}
  D₂::Form2{X}
  F::Form3{X}

  F == c₁(D₁)
  F == c₂(D₂)
end

Test1 = SummationDecapode(parse_decapode(DecaTest1))


Test1Res = rewrite_decapode(Test1)

# Test with multiple rewrites
# TODO: Support multiple different rewrites
DecaTest2 = quote
  C₁::Form0{X}
  C₂::Form0{X}
  D₁::Form1{X}
  D₂::Form1{X}
  F::Form1{X}

  D₁ == d₀(C₁)
  D₁ == d₀(C₂)
  F == c₁(D₁)
  F == c₂(D₂)

end

Test2 = SummationDecapode(parse_decapode(DecaTest2))

Test2Res = rewrite_decapode(Test2)

# Test to ensure isolation of rewrite from unsupported features
# TODO: Will be test for op1, op2,sum combined rewrite later
DecaTest3 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  E::Form2{X}
  F::Form3{X}
  G::Form4{X}

  G == ∧(A, B)
  G == k(C)
  G == t(D)
  G == F + E
end

Test3 = SummationDecapode(parse_decapode(DecaTest3))
Test3Res = rewrite_decapode(preTest3Res)

# Need to find a way to post-process the temps back off
post_match_3 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 2
  type = [:Form4, :Form4]
  name = [Symbol("••1"), :temp]

  Op1 = 1
  src = [2]
  tgt = [1]
  op1 = [:temp]
end

post_I_3 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 2
  type = [:Form4, :Form4]
  name = [Symbol("••1"), :temp]
end

post_sub_3 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 1
  type = [:Form4]
  name = [Symbol("••1")]
end

L = ACSetTransformation(post_I_3, post_match_3; Var = 1:2)
R = ACSetTransformation(post_I_3, post_sub_3; Var = [1, 1])

rule = Rule(L, R)

m = ACSetTransformation(post_match_3, Test3Res, Var = [5, 1], Op1 = [1])
postTest3Res = rewrite_match(rule, m)

# Test to ensure that ops from the same source are all preserved
DecaTest4 = quote
  C::Form9{X}
  D::Form4{X}

  D == k(C)
  D == t(C)
  D == p(C)
end

Test4 = SummationDecapode(parse_decapode(DecaTest4))
Test4Res = rewrite_decapode(Test4)

# Test that larger nary rewrites function properly
DecaTest5 = quote
  A::Form0{X}
  B::Form1{X}
  C::Form2{X}
  D::Form3{X}
  E::Form4{X}
  F::Form5{X}
  G::Form6{X}

  G == f(F)
  G == e(E)
  G == d(D)
  G == c(C)
  G == b(B)
  G == a(A)
end

Test5 = SummationDecapode(parse_decapode(DecaTest5))
# TODO: This rewrite takes a significant amount of time to complete
Test5Res = rewrite_decapode(Test5)

DecaTest6 = quote
  A::Form0{X}
  B::Form1{X}
  C::Form2{X}
  D::Form3{X}
  E::Form4{X}
  F::Form5{X}

  F == k(A)
  F == t(E)
  E == p(B, C)
  E == q(B, D)
end

Test6 = SummationDecapode(parse_decapode(DecaTest6))

Test6 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 6
  type = [:Form0, :Form1, :Form2, :Form3, :Form4, :Form5]
  name = [:A, :B, :C, :D, :E, :F]

  Op1 = 2
  src = [1, 5]
  tgt = [6, 6]
  op1 = [:k, :t]

  Op2 = 2
  proj1 = [2, 2]
  proj2 = [3, 4]
  res = [5, 5]
  op2 = [:p, :q]
end
#########################
# May be used later on to try out composing rewrite rules

# Source graph
G = @acset Graph begin V = 5; E = 4; src = [1, 2, 3, 5]; tgt = [3, 3, 4, 4]end

# Graphs represenative of what the decapode
# rewrite should be accomplishing

# To match for
Match = @acset Graph begin V = 3; E = 2; src = [1, 2]; tgt = [3, 3] end

# Change into
Sub = @acset Graph begin V = 6; E = 5; src = [1, 2, 4, 5, 6]; tgt = [5, 6, 3, 4, 4] end

# Preserved by rewrite
I = Graph(3)

L = CSetTransformation(I, Match, V = [1,2,3])
R = CSetTransformation(I, Sub, V = [1,2,3])

rule = Rule(oplus(L, L), oplus(R, R))

m = CSetTransformation(oplus(Match, Match), G, V=[1,2,3,3,5,4], E=[1,2,3,4])
H = rewrite_match(rule, m)

#H = rewrite(rule, G)

function makePerfectBinaryTree(h)
  Tree = Graph(2^h - 1)
  interiorNodes = 2^(h-1)-1
  add_edges!(Tree, map(x->2*x, 1:interiorNodes), 1:interiorNodes)
  add_edges!(Tree, map(x->2*x+1, 1:interiorNodes), 1:interiorNodes)
  return Tree
end

G′ = makePerfectBinaryTree(3)
m = CSetTransformation(Match, G′, V=[2,3,1], E=[1,4])
H′ = rewrite_match(rule, m)

m = CSetTransformation(Match, H′, V=[9,10,2], E=[7,9])
H′′ = rewrite_match(rule, m)

m = CSetTransformation(Match, H′′, V=[12,13,7], E=[11,12])
H′′′ = rewrite_match(rule, m)

function rewriteNAryGraph(n)
  L = @acset Graph begin
    V = n+1
    E = n
    src = 1:n
    tgt = n+1
  end

  R = Graph()
  copy_parts!(R, L)
  m = add_vertices!(R, n)
  add_edges!(R, m, 1:n)
  rightmost_vertex = add_vertex!(R)
  add_edge!(R, n+1, rightmost_vertex)

  I = Graph(n+1)

  L′ = CSetTransformation(I, L, V = 1:n+1)
  R′ = CSetTransformation(I, R, V = vcat(1:n, rightmost_vertex))

  rule = Rule(L′, R′)

  H = rewrite(rule, L)
end

