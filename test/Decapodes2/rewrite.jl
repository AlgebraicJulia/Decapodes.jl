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


#Source graph
G = @acset Graph begin
    V = 4
    E = 3
    src = [1, 2, 3]
    tgt = [3, 3, 4]
end

# Graphs represenative of what the decapode
# rewrite should be accomplishing

# To match for
L = @acset Graph begin
    V = 3
    E = 2
    src = [1, 2]
    tgt = [3, 3]
end

# Change into
R = @acset Graph begin
    V = 6
    E = 5
    src = [1, 2, 4, 5, 6]
    tgt = [5, 6, 3, 4, 4]
end

# Preserved by rewrite
I = @acset Graph begin
    V = 3
end

#rule = Rule(hom(L,L), hom(L,R))
#rule = Rule(hom(I,L), hom(I,I))

rule = Rule(hom(I,L), hom(I,R))

#Works for trivial pattern matching 
H = rewrite(rule, L)
#But not for pattern-finding
H₂ = rewrite(rule, G)


# Example from the rewrite documentation
G = @acset Graph begin
    V=3; E=3;
    src=[1,2,2];
    tgt=[2,3,3]
end

L = @acset Graph begin V=2; E=2; src=1; tgt=2 end # matched pattern
I = @acset Graph begin V=2; E=1; src=1; tgt=2 end # interface: non-deleted subset of L
R = @acset Graph begin V=1; E=1; src=1; tgt=1 end # Replacement pattern
rule = Rule(hom(I,L), hom(I,R))
H = rewrite(rule, G)