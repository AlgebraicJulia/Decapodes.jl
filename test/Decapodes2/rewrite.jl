using Decapodes
using Catlab.Graphics
using Catlab, Catlab.Graphs, Catlab.Graphics, Catlab.CategoricalAlgebra
using Catlab.Theories, Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite
const hom = AlgebraicRewriting.homomorphism

draw(g; kw...) = to_graphviz(g; node_labels=true, edge_labels=true, kw...)
draw(f::ACSetTransformation; kw...) =
  to_graphviz(f; node_labels=true, edge_labels=true, draw_codom=false, kw...)


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
  F == k(Z)
end

G = SummationDecapode(parse_decapode(DecaSource))

DecaI = quote
  C₁::Form0{X}
  C₂::Form0{X}
  Z::Form1{X}
end

I = SummationDecapode(parse_decapode(DecaI))

OpenSummationDecapodeOb, OpenSummationDecapode = OpenACSetTypes(SummationDecapode, :Var)

OpenSub = Open(Sub, [:C₁, :C₂, :Z])
OpenMatch = Open(Match, [:C₁, :C₂, :Z])
OpenG = Open(G, [:C₁, :C₂, :Z])
OpenI = Open(I, [:C₁, :C₂, :Z])

matchTrans = ACSetTransformation(I, Match, Var = [1,2,3]);
subTrans = ACSetTransformation(I, Sub, Var = [1,2,3]);

"""
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
"""

L_ = OpenSummationDecapodeOb.body

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
