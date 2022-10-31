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

# This could be made more efficient
# Largest nary stored first
function get_rewrite_count(deca_source)
  # Stores the nary of all possible rewrites locations
  nary_counts = []
  for var in parts(deca_source, :Var)
    op1Count = length(incident(deca_source, var, :tgt))
    # op2Count = length(incident(deca_source, var, :res))
    # sumCount = length(incident(deca_source, var, :sum))

    tot = op1Count # + op2Count + sumCount
    if(tot >= 2)
      push!(nary_counts, tot)
    end
  end

  return sort(nary_counts, rev = true)
end



function rewrite_decapode(deca_source)
  # Just for now, I'll get it working with op1 only
  # Considering op2 and summations will make this significantly more difficult

  nary_counts = get_rewrite_count(deca_source)
  deca_rewrite = deca_source

  for count in nary_counts
    nary_of_rewrite = first(nary_counts)
    num_nodes_match = nary_of_rewrite + 1
    result_index = num_nodes_match
    sum_index = 2 * result_index

    variable_types = map(n -> Var(Symbol("t",string(n))), 1:num_nodes_match)
    variable_var = map(n -> Var(Symbol("v",string(n))), 1:num_nodes_match)
    variable_op1 = map(n -> Var(Symbol("op1_",string(n))), 1:nary_of_rewrite)

    # Cannot work with undefined types
    Match_1 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
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

    I_1 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
      Var = num_nodes_match
      type = variable_types
      name = variable_var
    end 

    Sub_1 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
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

    L = hom(I_1, Match_1; monic = true, bindvars = true);
    R = hom(I_1, Sub_1; monic = true, bindvars = true);
    
    rule = Rule(L, R)

    deca_rewrite = rewrite(rule, deca_rewrite)
  end

  return deca_rewrite
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

Test1 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
  Var = 3
  type = [:Form1, :Form2, :Form3]
  name = [:D₁, :D₂, :F]

  Op1 = 2
  src = [1, 2]
  tgt = [3, 3]
  op1 = [:c₁, :c₂]
end

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
Test2 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
  Var = 5
  type = [:Form0, :Form0, :Form1, :Form1,:Form1,]
  name = [:C₁, :C₂, :D₁, :D₂, :F]

  Op1 = 4
  src = [1, 2, 3, 4]
  tgt = [3, 3, 5, 5]
  op1 = [:d₀, :d₀, :c₁, :c₂]
end

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

Test3 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
  Var = 7
  type = [:Form0, :Form0, :Form0, :Form0, :Form2, :Form3, :Form4]
  name = [:A, :B, :C, :D, :E, :F, :G]

  Op1 = 2
  src = [3, 4]
  tgt = [7, 7]
  op1 = [:k, :t]

  Op2 = 1
  proj1 = [1]
  proj2 = [2]
  res = [7]
  op2 = [:∧]

  Σ = 1
  sum = [7]

  Summand = 2
  summand = [6, 5]
  summation = [1, 1]
end

Test3Res = rewrite_decapode(Test3)

# Test to ensure that ops from the same source are all preserved
DecaTest4 = quote
  C::Form9{X}
  D::Form4{X}

  D == k(C)
  D == t(C)
  D == p(C)
end

Test4 = SummationDecapode(parse_decapode(DecaTest4))

Test4 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
  Var = 2
  type = [:Form9, :Form4]
  name = [:C, :D]

  Op1 = 3
  src = [1, 1, 1]
  tgt = [2, 2, 2]
  op1 = [:k, :t, :p]
end

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

Test5 = @acset SummationDecapode{Any, Any, Union{Var, Symbol}} begin 
  Var = 7
  type = [:Form0, :Form1, :Form2, :Form3, :Form4, :Form5, :Form6]
  name = [:A, :B, :C, :D, :E, :F, :G]

  Op1 = 6
  src = [6, 5, 4, 3, 2, 1]
  tgt = [7, 7, 7, 7, 7, 7]
  op1 = [:f, :e, :d, :c, :b, :a]
end

# TODO: This rewrite takes a significant amount of time to complete
Test5Res = rewrite_decapode(Test5)

#########################
# May be used later on to try out composing rewrite rules

# Source graph
G = @acset Graph begin V = 4; E = 3; src = [1, 2, 3]; tgt = [3, 3, 4]end

# Graphs represenative of what the decapode
# rewrite should be accomplishing

# To match for
Match = @acset Graph begin V = 3; E = 2; src = [1, 2]; tgt = [3, 3] end

# Change into
Sub = @acset Graph begin V = 6; E = 5; src = [1, 2, 4, 5, 6]; tgt = [5, 6, 3, 4, 4]
end

# Preserved by rewrite
I = Graph(3)

L = CSetTransformation(I, Match, V = [1,2,3])
R = CSetTransformation(I, Sub, V = [1,2,3])
rule = Rule(L, R)

#m = CSetTransformation(Match, G, V=[1,2,3], E=[1,2])
#H = rewrite_match(rule, m)

H = rewrite(rule, G)

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