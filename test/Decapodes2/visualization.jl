using Decapodes
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs.PropertyGraphs

DecaTest = quote
    A::Form0{X}
    B::Form0{X}
    C::Form0{X}
    D::Form0{X}
  
    D == k(A + B) + p(C + B)
  end
  
Test1 = SummationDecapode(parse_decapode(DecaTest))
to_graphviz(Test1)

DecaTest2 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  Ḋ::Form0{X}

  Ḋ == k(A, B) + p(C, B)
  ∂ₜ(D) == Ḋ 
end

Test2 = SummationDecapode(parse_decapode(DecaTest2))
to_graphviz(Test2)
to_graphviz(Test2, isDirected = false)