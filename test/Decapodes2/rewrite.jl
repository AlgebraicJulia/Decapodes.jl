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
tomatch = @acset Graph begin V = 3; E = 2; src = [1, 2]; tgt = [3, 3] end

# Change into
tosub = @acset Graph begin V = 6; E = 5; src = [1, 2, 4, 5, 6]; tgt = [5, 6, 3, 4, 4]
end

# Preserved by rewrite
I = Graph(3)

L = CSetTransformation(I, tomatch, V = [1,2,3])
R = CSetTransformation(I, tosub, V = [1,2,3])
rule = Rule(L, R)

#m = CSetTransformation(tomatch, G, V=[1,2,3], E=[1,2])
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

#H′ = rewrite(rule, G′)

ruleSeq = RuleSchedule(rule)
seq = WhileSchedule(ruleSeq)

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
tomatch = @acset Graph begin V = 3; E = 2; src = [1, 2]; tgt = [3, 3] end
tosub = @acset Graph begin V = 6; E = 5; src = [1, 2, 4, 5, 6]; tgt = [5, 6, 3, 4, 4] end

I = Graph(3)
id_1 = id(Graph(1));

# Create the open versions of each graph
openG = OpenGraph(G, FinFunction([1], 4), FinFunction([2], 4), FinFunction([3], 4))
openMatch = OpenGraph(tomatch, FinFunction([1], 3), FinFunction([2], 3), FinFunction([3], 3))
openSub = OpenGraph(tosub, FinFunction([1], 6), FinFunction([2], 6), FinFunction([3], 6))

openI = OpenGraph(I, FinFunction([1], 3), FinFunction([2], 3), FinFunction([3], 3))


# Create the equivalances between the matching and substitute graph
matchTrans = ACSetTransformation(I, tomatch, V = [1,2,3]);
subTrans = ACSetTransformation(I, tosub, V = [1,2,3]);

# Extend the transformation to the entire multicospan, including apex and feet
L = StructuredMultiCospanHom(openI, openMatch, ACSetTransformation[matchTrans, id_1, id_1, id_1])
R = StructuredMultiCospanHom(openI, openSub, ACSetTransformation[subTrans, id_1, id_1, id_1])

# Declare the rewrite rule using both transformations
rule = openrule(Span(L, R))

# Declare the matching transformation, can we let the program
# match the pattern automatically?
findTrans = ACSetTransformation(tomatch, G, V=[1,2,3], E=[1,2]);
m = StructuredMultiCospanHom(openMatch, openG, ACSetTransformation[findTrans, id_1, id_1, id_1])

#Rewrite and get the apex, which is the result
Hmcs = open_rewrite_match(rule, m)
H = apex(Hmcs)

#########################

draw(d) = to_graphviz(d)

DecaSub = quote
  C₁::Form0{X}
  C₂::Form0{X}
  R::Form1{X}

  # Fick's first law
  R ==  k(d₀(C₁) + d₀(C₂))
end

subExpr = parse_decapode(DecaSub)
Sub = SummationDecapode(subExpr)

DecaMatch = quote
  C₁::Form0{X}
  C₂::Form0{X}
  R::Form1{X}

  # Fick's first law
  R == d₀(C₁)
  R == d₀(C₂)
end

matchExpr = parse_decapode(DecaMatch)
Match = SummationDecapode(matchExpr)

OpenSummationDecapodeOb, OpenSummationDecapode = OpenACSetTypes(SummationDecapode, :Var)

OpenSub = Open(Sub, [:C₁, :C₂, :R])
OpenMatch = Open(Match, [:C₁, :C₂, :R])
