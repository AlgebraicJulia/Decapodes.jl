using Test
using Decapodes
using Catlab.Graphics
using Catlab.Graphics.Graphviz
using Catlab.Graphs.PropertyGraphs

DecaTest = quote
    (A, B, C, D)::Form0{X}
  
    D == k(A + B) + p(C + B)
  end
  
Test1 = SummationDecapode(parse_decapode(DecaTest))
t1 = to_graphviz(Test1)
@test Graphviz.filter_statements(t1, Graphviz.Edge, :label) == ["k",
  "p",
  "+",
  "+",
  "+"]
@test Graphviz.filter_statements(t1, Graphviz.Node, :label) == [
  "A:Ω₀",
  "B:Ω₀",
  "C:Ω₀",
  "D:Ω₀",
  "•1:Ω•",
  "sum_1:Ω•",
  "•2:Ω•",
  "sum_2:Ω•",
  "",
  "",
  "Σ1",
  "Σ2",
  "Σ3" ]
@test Graphviz.filter_statements(t1, Graphviz.Node, :shape) == [
  "none",
  "none",
  "circle",
  "circle",
  "circle"]

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
t2 = to_graphviz(Test2)
@test Graphviz.filter_statements(t2, Graphviz.Edge, :label) == ["∂ₜ", "π₁", "π₂", "k", "π₁", "π₂", "p", "+"]
@test Graphviz.filter_statements(t2, Graphviz.Node, :label) == ["A:Ω₀", "B:Ω₀", "C:Ω₀", "D:Ω₀", "Ḋ:Ω•", "•1:Ω•", "•2:Ω•", "", "", "Ω₀×Ω₀", "Ω₀×Ω₀", "Σ1"]
@test Graphviz.filter_statements(t2, Graphviz.Node, :shape) == ["none", "none", "rectangle", "rectangle", "circle"]

t2_undirected = to_graphviz(Test2, directed = false)
@test Graphviz.filter_statements(t2_undirected, Graphviz.Edge, :label) == ["∂ₜ", "π₁", "π₂", "k", "π₁", "π₂", "p", "+"]
@test Graphviz.filter_statements(t2_undirected, Graphviz.Node, :label) == ["A:Ω₀", "B:Ω₀", "C:Ω₀", "D:Ω₀", "Ḋ:Ω•", "•1:Ω•", "•2:Ω•", "Ω₀×Ω₀", "Ω₀×Ω₀", "Σ1"]
@test Graphviz.filter_statements(t2_undirected, Graphviz.Node, :shape) == ["rectangle", "rectangle", "circle"]
