using Decapodes
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs.PropertyGraphs

draw(g; kw...) = to_graphviz(g; node_labels=true, edge_labels=true, kw...)

DecaTest = quote
    (A, B, C, D)::Form0{X}
  
    D == k(A + B) + p(C + B)
  end
  
Test1 = SummationDecapode(parse_decapode(DecaTest))
