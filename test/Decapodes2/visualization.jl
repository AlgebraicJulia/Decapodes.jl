using Decapodes
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs.PropertyGraphs

draw(g; kw...) = to_graphviz(g; node_labels=true, edge_labels=true, kw...)

DecaTest = quote
    A::Form0{X}
    B::Form0{X}
    C::Form0{X}
    D::Form0{X}
  
    D == k(A + B) + p(C + B)
  end
  
Test1 = SummationDecapode(parse_decapode(DecaTest))