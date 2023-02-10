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
@test Graphviz.filter_statements(t1, Graphviz.Edge, :label) == ["k", "p", "+", "+", "+"]
@test Graphviz.filter_statements(t1, Graphviz.Node, :label) == [ "A:Ω₀", "B:Ω₀", "C:Ω₀", "D:Ω₀", "•1:Ω•", "sum_1:Ω•", "•2:Ω•", "sum_2:Ω•", "", "", "Σ1", "Σ2", "Σ3" ]
@test Graphviz.filter_statements(t1, Graphviz.Node, :shape) == [ "none", "none", "circle", "circle", "circle"]
# Test that the default attributes are correct.
@test t1.graph_attrs == Dict(:rankdir => "TB")
@test t1.node_attrs == Dict(:height => "0.05", :margin => "0", :shape => "oval", :width => "0.05")

t1_attributes = to_graphviz(Test1, edge_attrs=Dict(:color => "cornflowerblue"), node_attrs=Dict(:shape => "egg"), graph_attrs=Dict(:rankdir => "LR"))
# Test that the default attributes are overwritten only where specified.
@test t1_attributes.edge_attrs == Dict(:arrowsize => "0.5", :color => "cornflowerblue")
@test t1_attributes.node_attrs == Dict(:height => "0.05", :margin => "0", :shape => "egg", :width => "0.05")
@test t1_attributes.graph_attrs == Dict(:rankdir => "LR")
# Test that the per-edge and per-node attributes are not overwritten.
@test Graphviz.filter_statements(t1_attributes, Graphviz.Edge, :label) == ["k", "p", "+", "+", "+"]
@test Graphviz.filter_statements(t1_attributes, Graphviz.Node, :label) == [ "A:Ω₀", "B:Ω₀", "C:Ω₀", "D:Ω₀", "•1:Ω•", "sum_1:Ω•", "•2:Ω•", "sum_2:Ω•", "", "", "Σ1", "Σ2", "Σ3" ]
@test Graphviz.filter_statements(t1_attributes, Graphviz.Node, :shape) == [ "none", "none", "circle", "circle", "circle"]

t1_no_default_changes = to_graphviz(Test1, node_attrs=Dict(:color => "red"), graph_attrs=Dict(:bgcolor => "fuchsia"))
# Test that the default attributes are not overwritten when not specified.
# (In this case, there should be no overwriting of default node shape, and graph rankdir.)
@test t1_no_default_changes.node_attrs == Dict(:color => "red", :height => "0.05", :margin => "0", :shape => "oval", :width => "0.05")
@test t1_no_default_changes.graph_attrs == Dict(:bgcolor => "fuchsia", :rankdir => "TB")

DecaTest2 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  Ḋ::Form0{X}

  Ḋ == k(A, B) + p(C, B) - 1
  ∂ₜ(D) == Ḋ 
end

Test2 = SummationDecapode(parse_decapode(DecaTest2))
t2 = to_graphviz(Test2)
@test Graphviz.filter_statements(t2, Graphviz.Edge, :label) == ["∂ₜ", "π₁", "π₂", "k", "π₁", "π₂", "p", "π₁", "π₂", "-", "+"]
@test Graphviz.filter_statements(t2, Graphviz.Node, :label) == ["A:Ω₀", "B:Ω₀", "C:Ω₀", "D:Ω₀", "Ḋ:Ω₀", "1:ΩL", "•1:Ω•", "•2:Ω•", "sum_1:Ω•", "", "", "Ω₀×Ω₀", "Ω₀×Ω₀", "Ω•×ΩL", "Σ1"]
@test Graphviz.filter_statements(t2, Graphviz.Node, :shape) == ["none", "none", "rectangle", "rectangle", "rectangle", "circle"]

t2_undirected = to_graphviz(Test2, directed = false)
@test Graphviz.filter_statements(t2_undirected, Graphviz.Edge, :label) == ["∂ₜ", "π₁", "π₂", "k", "π₁", "π₂", "p", "π₁", "π₂", "-", "+"]
@test Graphviz.filter_statements(t2_undirected, Graphviz.Node, :label) == ["A:Ω₀", "B:Ω₀", "C:Ω₀", "D:Ω₀", "Ḋ:Ω₀", "1:ΩL", "•1:Ω•", "•2:Ω•", "sum_1:Ω•", "Ω₀×Ω₀", "Ω₀×Ω₀", "Ω•×ΩL", "Σ1"]
@test Graphviz.filter_statements(t2_undirected, Graphviz.Node, :shape) == ["rectangle", "rectangle", "rectangle", "circle"]
