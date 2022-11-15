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

function get_valid_op1s(deca_source, varID)
  # skip_ops = Set([:∂ₜ])
  indices = incident(deca_source, varID, :tgt)
  return filter!(x -> deca_source[x, :op1] != :∂ₜ, indices)
end

function get_target_indices(deca_source)
  targetVars = []
  for var in parts(deca_source, :Var)
    op1Count = length(get_valid_op1s(deca_source, var))
    op2Count = length(incident(deca_source, var, :res))
    sumCount = length(incident(deca_source, var, :sum))

    tot = op1Count + op2Count + sumCount
    if(tot >= 2)
      append!(targetVars, var)
    end
  end

  return targetVars
end

function get_preprocess_indices(deca_source)
  targetOp2 = []
  targetSum = []

  targetVars = get_target_indices(deca_source)

  for var in targetVars
    append!(targetOp2, incident(deca_source, var, :res))
    append!(targetSum, incident(deca_source, var, :sum))
  end

  return targetOp2, targetSum
end

function preprocess_rewrite(deca_source)
  targetOp2, targetSum = get_preprocess_indices(deca_source)

  # If we don't need to preprocess then don't
  if(length(targetOp2) == 0 && length(targetSum) == 0)
    return deca_source
  end

  LHS = []
  RHS = []

  SuperMatch = []
  SuperVarMap = Vector{Int64}()
  SuperOp2Map = Vector{Int64}()
  SuperSigmaMap = Vector{Int64}()
  SuperSummandMap = Vector{Int64}()

  serial = 0
  # Process all of the target rewrites for op2
  for opID in targetOp2

    vars = [deca_source[opID, :proj1], deca_source[opID, :proj2], deca_source[opID, :res]]
    types = deca_source[vars, :type]
    names = deca_source[vars, :name]
    op2name = deca_source[opID, :op2]

    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 3
      type = types
      name = names

      Op2 = 1
      proj1 = [1]
      proj2 = [2]
      res = [3]
      op2 = op2name
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 3
      type = types
      name = names
    end

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 4
      type = vcat(types, types[end])
      name = vcat(names, Symbol("Temp_", serial))

      Op1 = 1
      src = [4]
      tgt = [3]
      op1 = [:temp]

      Op2 = 1
      proj1 = [1]
      proj2 = [2]
      res = [4]
      op2 = op2name
    end

    serial += 1

    L = ACSetTransformation(I, Match, Var = 1:3)
    R = ACSetTransformation(I, Sub, Var = 1:3)

    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    push!(SuperOp2Map, opID)
  end
    
  # Process all of the target rewrites for sums
  for sumID in targetSum
    summandIDs = incident(deca_source, sumID, :summation)
    vars = vcat(deca_source[summandIDs, :summand], deca_source[sumID, :sum])
    types = deca_source[vars, :type]
    names = deca_source[vars, :name]
    
    rewrite_size = length(vars)
    nary = rewrite_size - 1

    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = rewrite_size
      type = types
      name = names

      Σ = 1
      sum = [rewrite_size]

      Summand = nary
      summand = 1:nary
      summation = fill(1, nary)
    end

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = rewrite_size
      type = types
      name = names
    end

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = rewrite_size + 1
      type = vcat(types, types[end])
      name = vcat(names, Symbol("Temp_", serial))

      Op1 = 1
      src = [rewrite_size + 1]
      tgt = [rewrite_size]
      op1 = [:temp]

      Σ = 1
      sum = [rewrite_size + 1]

      Summand = nary
      summand = 1:nary
      summation = fill(1, nary)
    end

    serial += 1

    L = ACSetTransformation(I, Match, Var = 1:rewrite_size)
    R = ACSetTransformation(I, Sub, Var = 1:rewrite_size)

    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    push!(SuperSigmaMap, sumID)
    append!(SuperSummandMap, summandIDs)
  end

  # Combine all rules in parallel and apply
  rule = Rule(oplus(LHS), oplus(RHS))
  m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op2 = SuperOp2Map, Σ = SuperSigmaMap, Summand = SuperSummandMap)

  rewrite_match(rule, m)
end

function rewrite_decapode(deca_source)
  # Just for now, I'll get it working with op1 only
  # Considering op2 and summations will make this significantly more difficult

  targetVars = get_target_indices(deca_source)

  if(length(targetVars) == 0)
    return deca_source
  end

  LHS = []
  RHS = []

  SuperMatch = []
  SuperVarMap = Vector{Int64}()
  SuperOp1Map = Vector{Int64}()

  varSerial = 0
  sumSerial = 0
  for varID in targetVars
    targetOp1 = get_valid_op1s(deca_source, varID)
    vars = vcat(deca_source[targetOp1, :src], varID)

    num_nodes_match = length(vars)
    nary_of_rewrite = num_nodes_match - 1

    result_index = num_nodes_match
    sum_index = 2 * result_index

    variable_types = deca_source[vars, :type]
    variable_var =  deca_source[vars, :name]
    variable_op1 = deca_source[targetOp1, :op1]

    # Cannot work with undefined types
    Match = @acset SummationDecapode{Any, Any, Symbol} begin 
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

    I = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = num_nodes_match
      type = variable_types
      name = variable_var
    end 

    Sub = @acset SummationDecapode{Any, Any, Symbol} begin 
      Var = 2 * num_nodes_match
      type = vcat(variable_types, variable_types)
      name = vcat(variable_var, map(x -> Symbol("••", varSerial + x), 1:nary_of_rewrite), [Symbol("••sum", sumSerial)])
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

    varSerial += nary_of_rewrite
    sumSerial += 1

    L = ACSetTransformation(I, Match, Var = 1:num_nodes_match);
    R = ACSetTransformation(I, Sub, Var = 1:num_nodes_match);
    
    push!(LHS, L)
    push!(RHS, R)
    push!(SuperMatch, Match)

    append!(SuperVarMap, vars)
    append!(SuperOp1Map, targetOp1)
  end

  rule = Rule(oplus(LHS), oplus(RHS))

  m = ACSetTransformation(oplus(SuperMatch), deca_source, Var = SuperVarMap, Op1 = SuperOp1Map)
  rewrite_match(rule, m)
end

# No valid rewrites, return original decapode
DecaTest0 = quote
  A::Form1{X}
  B::Form2{X}
  C::Form3{X}
  
  D::Form1{X}
  E::Form2{X}
  F::Form3{X}

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
Test0Res = rewrite_decapode(preprocess_rewrite(Test0))
@test Test0 == Test0Res

# Trivial rewrite test
# Make sure var, op names and forms are preserved
DecaTest1 = quote
  D₁::Form1{X}
  D₂::Form2{X}
  F::Form3{X}

  F == c₁(D₁)
  F == c₂(D₂)
end

Test1 = SummationDecapode(parse_decapode(DecaTest1))
Test1Res = rewrite_decapode(preprocess_rewrite(Test1))

Test1Expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 6
  type = [:Form1, :Form2, :Form3, :Form1, :Form2, :Form3]
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
Test2Res = rewrite_decapode(preprocess_rewrite(Test2))

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
  E::Form2{X}
  F::Form3{X}
  G::Form4{X}

  G == ∧(A, B)
  G == k(C)
  G == t(D)
  G == F + E
end

Test3 = SummationDecapode(parse_decapode(DecaTest3))
Test3Res = rewrite_decapode(preprocess_rewrite(Test3))

Test3Expected = @acset SummationDecapode{Any, Any, Symbol} begin
  Var = 14

  type = Any[:Form4, :Form4, :Form0, :Form0, :Form4, :Form4, 
  :Form4, :Form0, :Form0, :Form4, :Form0, :Form0, :Form3, :Form2]

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
  C::Form9{X}
  D::Form4{X}

  D == k(C)
  D == t(C)
  D == p(C)
end

Test4 = SummationDecapode(parse_decapode(DecaTest4))
Test4Res = rewrite_decapode(preprocess_rewrite(Test4))

Test4Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 6
  type = Any[:Form9, :Form4, :Form9, :Form9, :Form9, :Form4] 
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
Test5Res = rewrite_decapode(preprocess_rewrite(Test5))

Test5Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 14
  type = Any[:Form5, :Form4, :Form3, :Form2, :Form1, :Form0, 
  :Form6, :Form5, :Form4, :Form3, :Form2, :Form1, :Form0, :Form6]
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
  D::Form3{X}
  E::Form4{X}
  F::Form5{X}

  F == k(A)
  F == t(E)
  E == p(B, C)
  E == q(B, D)
end

Test6 = SummationDecapode(parse_decapode(DecaTest6))
Test6Res = rewrite_decapode(preprocess_rewrite(Test6))

Test6Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 14
  type = Any[:Form4, :Form4, :Form4, :Form4, :Form4, :Form4, 
  :Form0, :Form5, :Form0, :Form4, :Form5, :Form1, :Form2, :Form3]
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
  D::Form3{X}
  E::Form4{X}
  F::Form5{X}

  F == A + B
  F == C + D + E
end

Test7 = SummationDecapode(parse_decapode(DecaTest7))
Test7Res = rewrite_decapode(preprocess_rewrite(Test7))

Test7Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 11
  type = Any[:Form5, :Form5, :Form5, :Form5, :Form5, :Form5, 
  :Form0, :Form1, :Form2, :Form3, :Form4]
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
DecaTest8 = quote
  D₁::Form1{X}
  D₂::Form2{X}
  F::Form3{X}

  ∂ₜ(D₁) == F
  F == c₂(D₂)
end

Test8 = SummationDecapode(parse_decapode(DecaTest8))
Test8Res = rewrite_decapode(preprocess_rewrite(Test8))

Test8Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 4
  type = Any[:Form1, :Form2, :Form3, :infer]
  name = [:D₁, :D₂, :F, :D₁̇]

  TVar = 1
  incl = [4]

  Op1 = 2
  src = [1, 2]
  tgt = [4, 3]
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
Test9Res = rewrite_decapode(preprocess_rewrite(Test9))=#

# Test that rewrites preverse TVars, ignore ∂ₜ
DecaTest10 = quote
  A::Form0{X}
  B::Form0{X}
  C::Form0{X}
  D::Form0{X}
  Ḣ::Form5{X}
  H::Form5{X}

  A == b(B)
  A == d(D)
  Ḣ == c(C)
  ∂ₜ(H) == Ḣ
end

Test10 = SummationDecapode(parse_decapode(DecaTest10))
Test10Res = rewrite_decapode(preprocess_rewrite(Test10))

Test10Expected = @acset SummationDecapode{Any, Any, Symbol}  begin
  Var = 9
  type = Any[:Form0, :Form0, :Form0, :Form0, :Form0, :Form0, 
  :Form0, :infer, :Form5]
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
BinTestRes = rewrite_decapode(preprocess_rewrite(BinTest))

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
#=
#########################
# May be used later on to try out composing rewrite rules

# Source graph
G = @acset Graph begin V = 5; E = 4; src = [1, 2, 3, 5]; tgt = [3, 3, 4, 4]end

# Graphs represenative of what the decapode
# rewrite should be accomplishing

# To match for
Match = @acset Graph begin V = 3; E = 2; src = [1, 2]; tgt = [3, 3] end

# Change into
Sub = @acset Graph begin V = 6; E = 5; src = [1, 2, 4, 5, 6]; tgt = [5, 6, 3, 4, 4] end

# Preserved by rewrite
I = Graph(3)

L = CSetTransformation(I, Match, V = [1,2,3])
R = CSetTransformation(I, Sub, V = [1,2,3])

rule = Rule(oplus(L, L), oplus(R, R))

m = CSetTransformation(oplus(Match, Match), G, V=[1,2,3,3,5,4], E=[1,2,3,4])
H = rewrite_match(rule, m)

#H = rewrite(rule, G)

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


# Need to find a way to post-process the temps back off
# Seems to work if the temp being removed has no edges into it

TrialDeca = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 4
  type = [:Form0, :Form1, :Form2, :Form3]
  name = [:A, :B, :C, :D]

  Op1 = 3
  src = [1, 2, 3]
  tgt = [2, 3, 4]
  op1 = [:a, :b, :c]
end

ResDeca = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 3
  type = [:Form0, :Form2, :Form3]
  name = [:A, :C, :D]

  Op1 = 2
  src = [1, 2]
  tgt = [2, 3]
  op1 = [:a, :c,]
end

post_match_3 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 2
  type = [:Form1, :Form2]
  name = [:B, :C]

  Op1 = 1
  src = [1]
  tgt = [2]
  op1 = [:b]
end

post_I_3 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 1
  type = [:Form2]
  name = [:C]
end

post_sub_3 = @acset SummationDecapode{Any, Any, Symbol} begin 
  Var = 1
  type = [:Form2]
  name = [:C]
end

L = hom(post_I_3, post_match_3)
R = hom(post_I_3, post_sub_3)

rule = Rule(L, R)

m = hom(post_match_3, TrialDeca; monic = true)
postWrite = rewrite_match(rule, m)

=#

# Meant for debugging purposes only, gives decapode acset structure
# to make copy decapodes for testing quickly
function generate_test_decapode(deca_source)
  println("@acset SummationDecapode{Any, Any, Symbol}  begin")
  if(nparts(deca_source, :Var) > 0)
    println("Var = ", nparts(deca_source, :Var))
    println("type = ", deca_source[:type])
    println("name = ", deca_source[:name])
    println("")
  end

  if(nparts(deca_source, :TVar) > 0)
    println("TVar = ", nparts(deca_source, :TVar))
    println("incl = ", deca_source[:incl])
    println("")
  end

  if(nparts(deca_source, :Op1) > 0)
    println("Op1 = ", nparts(deca_source, :Op1))
    println("src = ", deca_source[:src])
    println("tgt = ", deca_source[:tgt])
    println("op1 = ", deca_source[:op1])
    println("")
  end

  if(nparts(deca_source, :Op2) > 0)
    println("Op2 = ", nparts(deca_source, :Op2))
    println("proj1 = ", deca_source[:proj1])
    println("proj2 = ", deca_source[:proj2])
    println("res = ", deca_source[:res])
    println("op2 = ", deca_source[:op2])
    println("")
  end

  if(nparts(deca_source, :Σ) > 0)
    println("Σ = ", nparts(deca_source, :Σ))
    println("sum = ", deca_source[:sum])
    println("")
  end

  if(nparts(deca_source, :Summand) > 0)
    println("Summand = ", nparts(deca_source, :Summand))
    println("summand = ", deca_source[:summand])
    println("summation = ", deca_source[:summation])
  end

  println("end")
end

