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
# Work on rewriting a decapode based only on Op1

function get_rewrite_indexes(decaSource)
  # Convention will be that the target vertex will be the final
  # value in the rewrite_indexes, all before are source vertices
  temp = map(target -> vcat(decaSource[:src][incident(decaSource, target, :tgt)], target), unique(decaSource[:tgt]))
  rewrite_indexes = []
  for indexes in temp
    # Check if the relation isn't a trivial substitution
    # Only one Op1 into target
    if(length(indexes) > 2)
      push!(rewrite_indexes, indexes)
    end
  end
  return rewrite_indexes
end

function rewrite_decapode(decaSource, decaToRewrite)
  # Just for now, I'll get it working with this example

  nary_of_rewrite = 2
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

  # This matching morphism needs to change if we are going to
  # rewrite on decapodes other than the original
  # m = ACSetTransformation(Match_1, decaToRewrite, Var=info, Op1=incident(decaSource, info[end], :tgt));
  # rewrite_match(rule, m)

  rewrite(rule, decaSource)
end

"""
Match_2 = @acset SummationDecapode{Any, Any, Symbol} begin 
Var = 3
type = [:Form1, :Form1, :Form1]
name = [:D₁, :D₂, :F]
Op1 = 2
src = [1, 2]
tgt = [3, 3]
op1 = [:c₁, :c₂]
end

I_2 = @acset SummationDecapode{Any, Any, Symbol} begin 
Var = 3
type = [:Form1, :Form1, :Form1]
name = [:D₁, :D₂, :F]
end 

Sub_2 = @acset SummationDecapode{Any, Any, Symbol} begin 
Var = 6
type = [:Form1, :Form1, :Form1, :infer, :infer, :infer]
name = [:D₁, :D₂, :F, :sum2, :Three, :Four]
Op1 = 3
src = [1,2,4]
tgt = [5,6,3]
op1 = [:c₁, :c₂, :k₂]
Σ = 1
sum = [4]
Summand = 2
summand = [5, 6]
summation = [1, 1]
end

L = ACSetTransformation(I_2, Match_2, Var = [1,2,3]);
R = ACSetTransformation(I_2, Sub_2, Var = [1,2,3]);

rule = Rule(L, R)

m = ACSetTransformation(Match_2, H, Var=[3,7,8], Op1=[4,5]);
H′′ = rewrite_match(rule, m)
"""

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
rewriteIndexesTest1 = get_rewrite_indexes(Test1)
@test rewriteIndexesTest1 == [[1,2,3]]

Test1Res = rewrite_decapode(Test1, Test1, rewriteIndexesTest1[1])

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

rewriteIndexesTest2 = get_rewrite_indexes(Test2)
@test rewriteIndexesTest2 == [[1,2,3], [3,4,5]]

# Should combine into one test once multiple rewrite working
Test2Res1 = rewrite_decapode(Test2, Test2, rewriteIndexesTest2[1])
Test2Res2 = rewrite_decapode(Test2, Test2, rewriteIndexesTest2[2])

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
rewriteIndexesTest3 = get_rewrite_indexes(Test3)
@test rewriteIndexesTest2 == [[3,4,7]]

TestRes3 = rewrite_decapode(Test3, Test3, rewriteIndexesTest3[1])

# Test to ensure that ops from the same source are all preserved
DecaTest4 = quote
  C::Form9{X}
  D::Form4{X}

  D == k(C)
  D == t(C)
  D == p(C)
end

Test4 = SummationDecapode(parse_decapode(DecaTest4))
rewriteIndexesTest4 = get_rewrite_indexes(Test4)
@test rewriteIndexesTest4 == [[1,1,1,2]]

TestRes4 = rewrite_decapode(Test4, Test4, rewriteIndexesTest4[1])

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
rewriteIndexesTest5 = get_rewrite_indexes(Test5)
@test rewriteIndexesTest5 == [[6, 5, 4, 3, 2, 1, 7]]

TestRes5 = rewrite_decapode(Test5, Test5, rewriteIndexesTest5[1])
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

""" This is rewritting using open types, 
don't know if I need this anymore.

#########################
const OpenGraphOb, OpenGraph = OpenCSetTypes(Graph, :V)

# Define the source graph, matching graph, substitute graph
G = @acset Graph begin V = 4; E = 3; src = [1, 2, 3]; tgt = [3, 3, 4]end
Match = @acset Graph begin V = 3; E = 2; src = [1, 2]; tgt = [3, 3] end
Sub = @acset Graph begin V = 6; E = 5; src = [1, 2, 4, 5, 6]; tgt = [5, 6, 3, 4, 4] end

I = Graph(3)
id_1 = id(Graph(1));

# Create the open versions of each graph
openG = OpenGraph(G, FinFunction([1], 4), FinFunction([2], 4), FinFunction([3], 4))
openMatch = OpenGraph(Match, FinFunction([1], 3), FinFunction([2], 3), FinFunction([3], 3))
openSub = OpenGraph(Sub, FinFunction([1], 6), FinFunction([2], 6), FinFunction([3], 6))

openI = OpenGraph(I, FinFunction([1], 3), FinFunction([2], 3), FinFunction([3], 3))


# Create the equivalances between the matching and substitute graph
matchTrans = ACSetTransformation(I, Match, V = [1,2,3]);
subTrans = ACSetTransformation(I, Sub, V = [1,2,3]);

# Extend the transformation to the entire multicospan, including apex and feet
L = StructuredMultiCospanHom(openI, openMatch, ACSetTransformation[matchTrans, id_1, id_1, id_1])
R = StructuredMultiCospanHom(openI, openSub, ACSetTransformation[subTrans, id_1, id_1, id_1])

# Declare the rewrite rule using both transformations
rule = openrule(Span(L, R))

# Declare the matching transformation, can we let the program
# match the pattern automatically?
findTrans = ACSetTransformation(Match, G, V=[1,2,3], E=[1,2]);
m = StructuredMultiCospanHom(openMatch, openG, ACSetTransformation[findTrans, id_1, id_1, id_1])

#Rewrite and get the apex, which is the result
Hmcs = open_rewrite_match(rule, m)
H = apex(Hmcs)

#########################
# Decapode rewriting using the OpenSummationDecapode 
# Not sure if it works as of now but leaving it here
OpenSummationDecapodeOb, OpenSummationDecapode = OpenACSetTypes(SummationDecapode, :Var)

DecaSub = quote
  C₁::Form0{X}
  C₂::Form0{X}
  Z::Form1{X}

  Z ==  k(d₀(C₁) + d₀(C₂))
end

Sub = SummationDecapode(parse_decapode(DecaSub))

DecaMatch = quote
  C₁::Form0{X}
  C₂::Form0{X}
  Z::Form1{X}

  Z == d₀(C₁)
  Z == d₀(C₂)
end

Match = SummationDecapode(parse_decapode(DecaMatch))

DecaSource = quote
  C₁::Form0{X}
  C₂::Form0{X}
  Z::Form1{X}
  F::Form1{X}

  Z == d₀(C₁)
  Z == d₀(C₂)
  F == c(Z)
end

G = SummationDecapode(parse_decapode(DecaSource))

DecaI = quote
  C₁::Form0{X}
  C₂::Form0{X}
  Z::Form1{X}
end

I = SummationDecapode(parse_decapode(DecaI))

OpenSub = Open(Sub, [:C₁, :C₂, :Z])
OpenMatch = Open(Match, [:C₁, :C₂, :Z])
OpenG = Open(G, [:C₁, :C₂, :Z])
OpenI = Open(I, [:C₁, :C₂, :Z])

matchTrans = ACSetTransformation(I, Match, Var = [1,2,3]);
subTrans = ACSetTransformation(I, Sub, Var = [1,2,3]);

DecaC1 = quote
  C₁::Form0{X}
end

DecaC2 = quote
  C₂::Form0{X}
end

DecaZ = quote
  Z::Form1{X}
end

C1 = SummationDecapode(parse_decapode(DecaC1))
C2 = SummationDecapode(parse_decapode(DecaC2))
Z = SummationDecapode(parse_decapode(DecaZ))

id_1 = id(Graph(1));

L_ = OpenSummationDecapodeOb{Symbol, Symbol, Symbol}.body

L = StructuredMultiCospanHom(OpenI, OpenMatch, ACSetTransformation[matchTrans, L_(id_1), L_(id_1), L_(id_1)])
R = StructuredMultiCospanHom(OpenI, OpenSub, ACSetTransformation[subTrans, id(C1), id(C2), id(Z)])

rule = openrule(Span(L, R))

findTrans = ACSetTransformation(Match, G, Var=[1,2,3], Op1=[1,2]);
m = StructuredMultiCospanHom(OpenMatch, OpenG, ACSetTransformation[findTrans, id(C1), id(C2), id(Z)])

# Crashing when it tries to invert a hom
# Seems to be caused because invert_hom seems to want a symbol passed to it

# Symbol for invert_hom seems to be retrieved from
# L_ = typeof(left(rule.data)).parameters[1]
# and then
# ob = L_.parameters[1]
# Where ob is the symbol to be used
# Instead of a symbol, like in the OpenGraph instance which gets a :V
# it gets the following, which is not a symbol,
# AnonACSet{TypeLevelBasicSchema{Symbol, Tuple{:Var}, Tuple{}, 
# Tuple{:Type, :Operator, :Name}, Tuple{(:type, :Var, :Type), (:name, 
# :Var, :Name)}}, Tuple{Any, Any, Symbol}, Catlab.LVectors.LVector{Int64,
# (:Var,)}, NamedTuple{(:type, :name), Tuple{Catlab.ColumnImplementations.
# DenseAttr{Any}, Catlab.ColumnImplementations.DenseAttr{Symbol}}}}

# Comparison
# OpenGraph (which works), typeof(left(rule.data))
# StructuredMultiCospanHom{Catlab.CategoricalAlgebra.StructuredCospans
# .FinSetDiscreteACSet{:V, Graph}}

#OpenPode (which doesn't work), typeof(left(rule.data))
# StructuredMultiCospanHom{Catlab.CategoricalAlgebra.StructuredCospans
#.DiscreteACSet{AnonACSet{TypeLevelBasicSchema{Symbol, Tuple{:Var}, 
# Tuple{}, Tuple{:Type, :Operator, :Name}, Tuple{(:type, :Var, :Type), 
# (:name, :Var, :Name)}}, Tuple{Any, Any, Symbol}, Catlab.LVectors
# .LVector{Int64, (:Var,)}, NamedTuple{(:type, :name), Tuple{Catlab
# .ColumnImplementations.DenseAttr{Any}, Catlab.ColumnImplementations
# .DenseAttr{Symbol}}}}, SummationDecapode{Any, Any, Symbol}}}

# Relevent code appears to be in AlgebraicRewriting>src>
# StructuredCospans>open_pushout_complement 
Hmcs = open_rewrite_match(rule, m)
H = apex(Hmcs)

"""