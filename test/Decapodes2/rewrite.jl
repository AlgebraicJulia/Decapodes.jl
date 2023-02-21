using Test

using Decapodes
using Catlab.CSetDataStructures
import Decapodes: average_rewrite

#= draw(g; kw...) = to_graphviz(g; node_labels=true, edge_labels=true, kw...)
draw(f::ACSetTransformation; kw...) =
  to_graphviz(f; node_labels=true, edge_labels=true, draw_codom=false, kw...) =#

#########################

# No valid rewrites, return original decapode
DecaTest0 = quote
  A::Form0{X}
  B::Form1{X}
  C::Form2{X}
  
  D::Form0{X}
  E::Form1{X}
  F::Form2{X}

  G::Form0{X}
  H::Form0{X}
  I::Form0{X}
  J::Form0{X}

  B == a(A)
  C == b(B)

  F == f(D, E)

  J == G + H + I
end

Test0 = SummationDecapode(parse_decapode(DecaTest0))
Test0Res = average_rewrite(Test0)
@test Test0 == Test0Res

# Trivial rewrite test
# Make sure var, op names and forms are preserved
DecaTest1 = quote
  D₁::Form0{X}
  D₂::Form1{X}
  F::Form2{X}

  F == c₁(D₁)
  F == c₂(D₂)
end

Test1 = SummationDecapode(parse_decapode(DecaTest1))
Test1Res = average_rewrite(Test1)

Test1Expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 6
  type = [:Form0, :Form1, :Form2, :Form0, :Form1, :Form2]
  name = [:D₁, :D₂, :F, Symbol("••1"), Symbol("••2"), Symbol("••sum0")]

  Op1 = 3
  src = [1, 2, 6]
  tgt = [4, 5, 3]
  op1 = [:c₁, :c₂, :avg2]

  Σ = 1
  sum = [6]

  Summand = 2
  summand = [4, 5]
  summation = [1, 1]
end

@test Test1Res == Test1Expected

# Test with multiple rewrites, op1 only
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
Test2Res = average_rewrite(Test2)

Test2Expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 11

  type = [:Form0, :Form0, :Form1, :Form0, :Form0, 
  :Form1, :Form1, :Form1, :Form1, :Form1, :Form1]

  name = [:C₁, :C₂, :D₁, Symbol("••1"), Symbol("••2"), Symbol("••sum0"), 
  :D₂, :F, Symbol("••3"), Symbol("••4"), Symbol("••sum1")]

  Op1 = 6
  src = [1, 2, 6, 3, 7, 11]
  tgt = [4, 5, 3, 9, 10, 8]
  op1 = [:d₀, :d₀, :avg2, :c₁, :c₂, :avg2]

  Σ = 2
  sum = [6, 11]

  Summand = 4
  summand = [4, 5, 9, 10]
  summation = [1, 1, 2, 2]
end

@test Test2Res == Test2Expected

# Test to ensure rewrites work for op1, op2, and sums
DecaTest3 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  E::Form1{X}
  F::Form2{X}
  G::Form2{X}

  G == ∧(A, B)
  G == k(C)
  G == t(D)
  G == F + E
end

Test3 = SummationDecapode(parse_decapode(DecaTest3))
Test3Res = average_rewrite(Test3)

Test3Expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 14

  type = Any[:Form2, :Form2, :Form0, :Form0, :Form2, :Form2, 
  :Form2, :Form0, :Form0, :Form2, :Form0, :Form0, :Form2, :Form1]

  name = [:Temp_0, :Temp_1, :C, :D, :G, Symbol("••1"), Symbol("••2"), 
  Symbol("••3"), Symbol("••4"), Symbol("••sum0"), :A, :B, :F, :E]

  Op1 = 5
  src = [1, 2, 3, 4, 10]
  tgt = [6, 7, 8, 9, 5]
  op1 = [:temp, :temp, :k, :t, :avg4]

  Op2 = 1
  proj1 = [11]
  proj2 = [12]
  res = [1]
  op2 = [:∧]

  Σ = 2
  sum = [10, 2]

  Summand = 6
  summand = [6, 7, 8, 9, 13, 14]
  summation = [1, 1, 1, 1, 2, 2]
end

@test Test3Res == Test3Expected

# Test to ensure that ops from the same source are all preserved
DecaTest4 = quote
  C::Form0{X}
  D::Form1{X}

  D == k(C)
  D == t(C)
  D == p(C)
end

Test4 = SummationDecapode(parse_decapode(DecaTest4))
Test4Res = average_rewrite(Test4)

Test4Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 6
  type = Any[:Form0, :Form1, :Form0, :Form0, :Form0, :Form1] 
  name = [:C, :D, Symbol("••1"), Symbol("••2"), Symbol("••3"), Symbol("••sum0")]

  Op1 = 4
  src = [1, 1, 1, 6]
  tgt = [3, 4, 5, 2]
  op1 = Any[:k, :t, :p, :avg3]

  Σ = 1
  sum = [6]

  Summand = 3
  summand = [3, 4, 5]
  summation = [1, 1, 1]
end

@test Test4Res == Test4Expected

# Test that larger nary rewrites function properly
DecaTest5 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form1{X}
  D::Form2{X}
  E::Form1{X}
  F::Form0{X}
  G::Form2{X}

  G == f(F)
  G == e(E)
  G == d(D)
  G == c(C)
  G == b(B)
  G == a(A)
end

Test5 = SummationDecapode(parse_decapode(DecaTest5))
Test5Res = average_rewrite(Test5)

Test5Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 14
  type = Any[:Form0, :Form1, :Form2, :Form1, :Form0, :Form0, 
  :Form2, :Form0, :Form1, :Form2, :Form1, :Form0, :Form0, :Form2]
  name = [:F, :E, :D, :C, :B, :A, :G, Symbol("••1"), Symbol("••2"), Symbol("••3"), Symbol("••4"), Symbol("••5"), Symbol("••6"), Symbol("••sum0")]

  Op1 = 7
  src = [1, 2, 3, 4, 5, 6, 14]
  tgt = [8, 9, 10, 11, 12, 13, 7]
  op1 = Any[:f, :e, :d, :c, :b, :a, :avg6]

  Σ = 1
  sum = [14]

  Summand = 6
  summand = [8, 9, 10, 11, 12, 13]
  summation = [1, 1, 1, 1, 1, 1]
end

@test Test5Res == Test5Expected

# Test multiple rewrites with op2
DecaTest6 = quote
  A::Form0{X}
  B::Form1{X}
  C::Form2{X}
  D::Form0{X}
  E::Form1{X}
  F::Form2{X}

  F == k(A)
  F == t(E)
  E == p(B, C)
  E == q(B, D)
end

Test6 = SummationDecapode(parse_decapode(DecaTest6))
Test6Res = average_rewrite(Test6)

Test6Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 14
  type = Any[:Form1, :Form1, :Form1, :Form1, :Form1, :Form1, 
  :Form0, :Form2, :Form0, :Form1, :Form2, :Form1, :Form2, :Form0]
  name = [:Temp_0, :Temp_1, :E, Symbol("••1"), Symbol("••2"), Symbol("••sum0"), :A, :F, Symbol("••3"), Symbol("••4"), Symbol("••sum1"), :B, :C, :D]

  Op1 = 6
  src = [1, 2, 6, 7, 3, 11]
  tgt = [4, 5, 3, 9, 10, 8]
  op1 = Any[:temp, :temp, :avg2, :k, :t, :avg2]

  Op2 = 2
  proj1 = [12, 12]
  proj2 = [13, 14]
  res = [1, 2]
  op2 = Any[:p, :q]

  Σ = 2
  sum = [6, 11]

  Summand = 4
  summand = [4, 5, 9, 10]
  summation = [1, 1, 2, 2]
end

@test Test6Res == Test6Expected

# Test multiple rewrites with sums
DecaTest7 = quote
  A::Form0{X}
  B::Form1{X}
  C::Form2{X}
  D::Form1{X}
  E::Form0{X}
  F::Form2{X}

  F == A + B
  F == C + D + E
end

Test7 = SummationDecapode(parse_decapode(DecaTest7))
Test7Res = average_rewrite(Test7)

Test7Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 11
  type = Any[:Form2, :Form2, :Form2, :Form2, :Form2, :Form2, 
  :Form0, :Form1, :Form2, :Form1, :Form0]
  name = [:Temp_0, :Temp_1, :F, Symbol("••1"), Symbol("••2"), Symbol("••sum0"), :A, :B, :C, :D, :E]

  Op1 = 3
  src = [1, 2, 6]
  tgt = [4, 5, 3]
  op1 = Any[:temp, :temp, :avg2]

  Σ = 3
  sum = [6, 1, 2]

  Summand = 7
  summand = [4, 5, 7, 8, 9, 10, 11]
  summation = [1, 1, 2, 2, 3, 3, 3]
end

@test Test7Res == Test7Expected

# Test that rewrite ignores forbidden ops, like ∂ₜ
# TODO: This test might break if TVar behavior changes
DecaTest8 = quote 
  D₁::Form1{X}
  D₂::Form2{X}
  F::Form0{X}

  ∂ₜ(D₁) == F
  F == c₂(D₂)
end

Test8 = SummationDecapode(parse_decapode(DecaTest8))
Test8Res = average_rewrite(Test8)

Test8Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 3
  type = Any[:Form1, :Form2, :Form0]
  name = [:D₁, :D₂, :D₁̇ ]

  TVar = 1
  incl = [3]

  Op1 = 2
  src = [1, 2]
  tgt = [3, 3]
  op1 = Any[:∂ₜ, :c₂]
end

@test Test8Res == Test8Expected

#=DecaTest9 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  E::Form2{X}
  F::Form3{X}
  Ḣ::Form5{X}
  H::Form5{X}

  Ḣ == k(B)
  Ḣ == p(D)
  ∂ₜ(H) == Ḣ
end

Test9 = SummationDecapode(parse_decapode(DecaTest9))
Test9Res = average_rewrite(Test9)=#

# Test that rewrites preverse TVars, ignore ∂ₜ
DecaTest10 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  Ḣ::Form2{X}
  H::Form2{X}

  A == b(B)
  A == d(D)
  Ḣ == c(C)
  ∂ₜ(H) == Ḣ
end

Test10 = SummationDecapode(parse_decapode(DecaTest10))
Test10Res = average_rewrite(Test10)

Test10Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 9
  type = Any[:Form0, :Form0, :Form0, :Form0, :Form0, :Form0, 
  :Form0, :Form2, :Form2]
  name = [:B, :D, :A, Symbol("••1"), Symbol("••2"), Symbol("••sum0"), :C, :Ḣ, :H]

  TVar = 1
  incl = [8]

  Op1 = 5
  src = [1, 2, 6, 7, 9]
  tgt = [4, 5, 3, 8, 8]
  op1 = Any[:b, :d, :avg2, :c, :∂ₜ]

  Σ = 1
  sum = [6]

  Summand = 2
  summand = [4, 5]
  summation = [1, 1]
end

@test Test10Res == Test10Expected

# Test for benchmarking large results, still gives correct output
function makePerfectBinaryDeca(h)
  num_nodes = 2^h - 1
  num_interior = 2^(h-1) - 1
  BinaryTest = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = num_nodes
    type = fill(:Form1, num_nodes)
    name = map(x -> Symbol("•$x"), 1:num_nodes)

    Op1 = 2 * num_interior
    src = vcat(map(x->2*x, 1:num_interior), map(x->2*x+1, 1:num_interior))
    tgt = vcat(1:num_interior, 1:num_interior)
    op1 = map(x -> Symbol("op$x"), 1:(2 * num_interior))

  end
end

h = 5
BinTest = makePerfectBinaryDeca(h)
BinTestRes = average_rewrite(BinTest)

BinTestExpected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 76
  type = Any[:Form1, :Form1, :Form1, :Form1, :Form1, :Form1, 
  :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, 
  :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1, :Form1]
  name = [Symbol("•2"), Symbol("•3"), Symbol("•1"), Symbol("••1"), Symbol("••2"), Symbol("••sum0"), Symbol("•4"), Symbol("•5"), Symbol("••3"), Symbol("••4"), Symbol("••sum1"), Symbol("•6"), Symbol("•7"), Symbol("••5"), Symbol("••6"), Symbol("••sum2"), Symbol("•8"), Symbol("•9"), Symbol("••7"), Symbol("••8"), Symbol("••sum3"), Symbol("•10"), Symbol("•11"), Symbol("••9"), Symbol("••10"), Symbol("••sum4"), Symbol("•12"), Symbol("•13"), Symbol("••11"), Symbol("••12"), Symbol("••sum5"), Symbol("•14"), Symbol("•15"), Symbol("••13"), 
  Symbol("••14"), Symbol("••sum6"), Symbol("•16"), Symbol("•17"), Symbol("••15"), Symbol("••16"), Symbol("••sum7"), Symbol("•18"), Symbol("•19"), Symbol("••17"), Symbol("••18"), Symbol("••sum8"), Symbol("•20"), Symbol("•21"), Symbol("••19"), Symbol("••20"), Symbol("••sum9"), Symbol("•22"), Symbol("•23"), Symbol("••21"), Symbol("••22"), Symbol("••sum10"), Symbol("•24"), Symbol("•25"), Symbol("••23"), Symbol("••24"), Symbol("••sum11"), Symbol("•26"), Symbol("•27"), Symbol("••25"), Symbol("••26"), Symbol("••sum12"), Symbol("•28"), Symbol("•29"), Symbol("••27"), Symbol("••28"), Symbol("••sum13"), Symbol("•30"), Symbol("•31"), Symbol("••29"), Symbol("••30"), Symbol("••sum14")]

  Op1 = 45
  src = [1, 2, 6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23, 26, 27, 28, 31, 32, 33, 36, 37, 38, 41, 42, 43, 46, 47, 48, 
  51, 52, 53, 56, 57, 58, 61, 62, 63, 66, 67, 68, 71, 72, 73, 76]
  tgt = [4, 5, 3, 9, 10, 1, 14, 15, 2, 19, 20, 7, 24, 25, 8, 
  29, 30, 12, 34, 35, 13, 39, 40, 17, 44, 45, 18, 49, 50, 22, 54, 55, 23, 59, 60, 27, 64, 65, 28, 69, 70, 32, 74, 75, 33]
  op1 = Any[:op1, :op16, :avg2, :op2, :op17, :avg2, :op3, :op18, :avg2, :op4, :op19, :avg2, :op5, :op20, :avg2, :op6, :op21, :avg2, :op7, :op22, :avg2, :op8, :op23, :avg2, :op9, :op24, :avg2, :op10, :op25, :avg2, :op11, :op26, :avg2, :op12, :op27, :avg2, :op13, :op28, :avg2, :op14, :op29, :avg2, 
  :op15, :op30, :avg2]

  Σ = 15
  sum = [6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76]

  Summand = 30
  summand = [4, 5, 9, 10, 14, 15, 19, 20, 24, 25, 29, 30, 34, 35, 39, 40, 44, 45, 49, 50, 54, 55, 59, 60, 64, 65, 69, 70, 74, 75]
  summation = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15]    
end

@test BinTestRes == BinTestExpected
