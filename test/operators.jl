@testset "Open Operators" begin
  laplace_de_rham_0 = @decapode begin
    (A, B)::Form0

    B == Δ₀(A)
  end
  test_compare = Decapodes.open_operators(laplace_de_rham_0)
  infer_types!(test_compare)
  resolve_overloads!(test_compare)

  Decapodes.open_operators!(laplace_de_rham_0)
  infer_types!(laplace_de_rham_0)
  resolve_overloads!(laplace_de_rham_0)


  @test test_compare == laplace_de_rham_0
  @test test_compare !== laplace_de_rham_0

end

@testset "Opening 1D Operators" begin
  # Test for 1D Interior Product (Dual1, Primal1 -> Dual0)
  interior_product_1 = @decapode begin
    A::Form1
    B::DualForm1
    C::DualForm0

    C == i₁(A, B)
  end
  Decapodes.open_operators!(interior_product_1; dimension = 1)
  infer_types!(interior_product_1, dim=1)
  resolve_overloads!(interior_product_1, dim=1)

  test_interior_product_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 5
    Op1 = 2
    Op2 = 1
    src = [3, 1]
    tgt = [2, 4]
    proj1 = [2]
    proj2 = [5]
    res = [1]
    op1 = [:⋆₀⁻¹, :⋆₁]
    op2 = [:∧₀₁]
    type = [:Form1, :Form0, :DualForm1, :DualForm0, :Form1]
    name = [Symbol("•1"), Symbol("•2"), :B, :C, :A]
  end
  @test interior_product_1 == test_interior_product_1

  # Test for 1D Lie Derivative (Dual0, Primal1 -> Dual0)
  lie_derivative_0 = @decapode begin
    A::DualForm0
    B::Form1
    C::DualForm0

    C == L₀(B, A)
  end
  Decapodes.open_operators!(lie_derivative_0; dimension = 1)
  infer_types!(lie_derivative_0, dim=1)
  resolve_overloads!(lie_derivative_0, dim=1)

  test_lie_derivative_0 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 6
    Op1 = 3
    Op2 = 1
    src = [1, 4, 2]
    tgt = [4, 3, 5]
    proj1 = [3]
    proj2 = [6]
    res = [2]
    op1 = [:dual_d₀, :⋆₀⁻¹, :⋆₁]
    op2 = [:∧₀₁]
    type = [:DualForm0, :Form1, :Form0, :DualForm1, :DualForm0, :Form1]
    name = [:A, Symbol("•1"), Symbol("•2"), Symbol("•1"), :C, :B]
  end
  @test lie_derivative_0 == test_lie_derivative_0

  # Test for 1D Lie Derivative (Dual1, Primal1 -> Dual1)
  lie_derivative_1 = @decapode begin
    A::DualForm1
    B::Form1
    C::DualForm1

    C == L₁(B, A)
  end
  Decapodes.open_operators!(lie_derivative_1; dimension = 1)
  infer_types!(lie_derivative_1, dim=1)
  resolve_overloads!(lie_derivative_1, dim=1)

  test_lie_derivative_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 6
    Op1 = 3
    Op2 = 1
    src = [5, 3, 1]
    tgt = [4, 2, 5]
    proj1 = [2]
    proj2 = [6]
    res = [1]
    op1 = [:dual_d₀, :⋆₀⁻¹, :⋆₁]
    op2 = [:∧₀₁]
    type = [:Form1, :Form0, :DualForm1, :DualForm1, :DualForm0, :Form1]
    name = [Symbol("•1"), Symbol("•2"), :A, :C, Symbol("•1"), :B]
  end
  @test lie_derivative_1 == test_lie_derivative_1

  # Test for Codifferential (Primal1 -> Primal0)
  codiff_1 = @decapode begin
    A::Form1
    B::Form0

    B == δ₁(A)
  end
  Decapodes.open_operators!(codiff_1, dimension = 1)
  infer_types!(codiff_1, dim=1)
  resolve_overloads!(codiff_1, dim=1)

  test_codiff_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 4
    Op1 = 3
    src = [4, 1, 3]
    tgt = [2, 3, 4]
    op1 = [:⋆₀⁻¹, :⋆₁, :dual_d₀]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    type = [:Form1, :Form0, :DualForm0, :DualForm1]
    name = [:A, :B, Symbol("•_1_1"), Symbol("•_1_2")]
  end
  @test codiff_1 == test_codiff_1
  
  # Test for Laplace de Rham (Primal0 -> Primal0)
  laplace_de_rham_0 = @decapode begin
    (A, B)::Form0

    B == Δ₀(A)
  end
  Decapodes.open_operators!(laplace_de_rham_0, dimension = 1)
  infer_types!(laplace_de_rham_0, dim=1)
  resolve_overloads!(laplace_de_rham_0, dim=1)

  test_laplace_de_rham_0 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 5
    Op1 = 4
    src = [5, 1, 2, 4]
    tgt = [3, 2, 4, 5]
    op1 = [:⋆₀⁻¹, :d₀, :⋆₁, :dual_d₀]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    type = [:Form0, :Form1, :Form0, :DualForm0, :DualForm1]
    name = [:A, Symbol("•_1_1"), :B, Symbol("•_1_1"), Symbol("•_1_2")]
  end
  @test laplace_de_rham_0 == test_laplace_de_rham_0

  # Test for Laplace de Rham (Primal1 -> Primal1)
  laplace_de_rham_1 = @decapode begin
    (A, B)::Form1

    B == Δ₁(A)
  end
  Decapodes.open_operators!(laplace_de_rham_1, dimension = 1)
  infer_types!(laplace_de_rham_1, dim=1)
  resolve_overloads!(laplace_de_rham_1, dim=1)

  test_laplace_de_rham_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 5
    type = [:Form1, :Form1, :Form0, :DualForm0, :DualForm1]    
    Op1 = 4
    src = [3, 5, 1, 4]
    tgt = [2, 3, 4, 5]
    op1 = [:d₀, :⋆₀⁻¹, :⋆₁, :dual_d₀]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    name = [:A, :B, Symbol("•_1_1"), Symbol("•_2_1"), Symbol("•_2_2")]
  end  
  @test laplace_de_rham_1 == test_laplace_de_rham_1
end

@testset "Opening 2D Operators" begin
  # Test for Interior Product (Dual1, Primal1 -> Dual0)
  interior_product_1 = @decapode begin
    A::Form1
    B::DualForm1
    C::DualForm0

    C == i₁(A, B)
  end
  Decapodes.open_operators!(interior_product_1)
  infer_types!(interior_product_1)
  resolve_overloads!(interior_product_1)

  test_interior_product_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 7
    Op1 = 2
    Op2 = 2
    src = [3, 1]
    tgt = [2, 7]
    proj1 = [6, 2]
    proj2 = [7, 5]
    res = [4, 1]
    op1 = [:⋆₁⁻¹, :⋆₂]
    op2 = [:*, :∧₁₁]
    type = [:Form2, :Form1, :DualForm1, :DualForm0, :Form1, :Literal, :DualForm0]
    name = [Symbol("•2"), Symbol("•3"), :B, :C, :A, Symbol("-1"), Symbol("•1")]
  end
  @test interior_product_1 == test_interior_product_1

  # Test for Interior Product (Dual2, Primal1 -> Dual1)
  interior_product_2 = @decapode begin
    A::Form1
    B::DualForm2
    C::DualForm1

    C == i₂(A, B)
  end
  Decapodes.open_operators!(interior_product_2)
  infer_types!(interior_product_2)
  resolve_overloads!(interior_product_2)

  test_interior_product_2 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 3
    name = [:A, :B, :C]
    type = [:Form1, :DualForm2, :DualForm1]
    Op1 = 0
    src = Int64[]
    tgt = Int64[]
    op1 = Any[]
    Op2 = 1
    proj1 = [1]
    proj2 = [2]
    res = [3]
    op2 = [:i₂]
  end
  @test interior_product_2 == test_interior_product_2

  # Test for Lie Derivative (Primal1, Dual0 -> Dual0)
  lie_derivative_0 = @decapode begin
    A::DualForm0
    B::Form1
    C::DualForm0

    C == L₀(B, A)
  end
  Decapodes.open_operators!(lie_derivative_0)
  infer_types!(lie_derivative_0)
  resolve_overloads!(lie_derivative_0)

  test_lie_derivative_0 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 8
    name = [:A, Symbol("•2"), Symbol("•3"), Symbol("•1"), :C, :B, Symbol("-1"), Symbol("•1")]
    type = [:DualForm0, :Form2, :Form1, :DualForm1, :DualForm0, :Form1, :Literal, :DualForm0]
    Op1 = 3
    src = [1, 4, 2]
    tgt = [4, 3, 8]
    op1 = [:dual_d₀, :⋆₁⁻¹, :⋆₂]
    Op2 = 2
    proj1 = [7, 3]
    proj2 = [8, 6]
    res = [5, 2]
    op2 = [:*, :∧₁₁]
  end
  @test lie_derivative_0 == test_lie_derivative_0

  # Test for Lie Derivative (Primal1, Dual1 -> Dual1)
  lie_derivative_1 = @decapode begin
    A::DualForm1
    B::Form1
    C::DualForm1

    C == L₁(B, A)
  end
  Decapodes.open_operators!(lie_derivative_1)
  infer_types!(lie_derivative_1)
  resolve_overloads!(lie_derivative_1)

  test_lie_derivative_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 11
    name = [Symbol("•2"), Symbol("•3"), Symbol("•3"), :C, Symbol("•1"), :A, Symbol("•2"), Symbol("•4"), :B, Symbol("-1"), Symbol("•1")]
    type = [:Form2, :DualForm1, :Form1, :DualForm1, :DualForm1, :DualForm1, :DualForm2, :DualForm0, :Form1, :Literal, :DualForm0]
    Op1 = 4
    src = [6, 8, 6, 1]
    tgt = [7, 2, 3, 11]
    op1 = [:dual_d₁, :dual_d₀, :⋆₁⁻¹, :⋆₂]
    Op2 = 3
    proj1 = [10, 9, 3]
    proj2 = [11, 7, 9]
    res = [8, 5, 1]
    op2 = [:*, :i₂, :∧₁₁]
    Σ = 1
    sum = [4]
    Summand = 2
    summand = [5, 2]
    summation = [1, 1]
  end
  @test lie_derivative_1 == test_lie_derivative_1

  # Test for Lie Derivative (Primal1, Dual2 -> Dual2)
  lie_derivative_2 = @decapode begin
    A::DualForm2
    B::Form1
    C::DualForm2

    C == L₂(B, A)
  end
  Decapodes.open_operators!(lie_derivative_2)
  infer_types!(lie_derivative_2)
  resolve_overloads!(lie_derivative_2)

  test_lie_derivative_2 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 4
    name = [:A, Symbol("•1"), :B, :C]
    type = [:DualForm2, :DualForm1, :Form1, :DualForm2]
    Op1 = 1
    src = [2]
    tgt = [4]
    op1 = [:dual_d₁]
    Op2 = 1
    proj1 = [3]
    proj2 = [1]
    res = [2]
    op2 = [:i₂]
  end
  @test lie_derivative_2 == test_lie_derivative_2

  # Test for Codifferential (Primal1 -> Primal0)
  codiff_1 = @decapode begin
    A::Form1
    B::Form0

    B == δ₁(A)
  end
  Decapodes.open_operators!(codiff_1)
  infer_types!(codiff_1)
  resolve_overloads!(codiff_1)

  test_codiff_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 4
    name = [:A, :B, Symbol("•_1_1"), Symbol("•_1_2")]
    type = [:Form1, :Form0, :DualForm1, :DualForm2]
    Op1 = 3
    src = [4, 1, 3]
    tgt = [2, 3, 4]
    op1 = [:⋆₀⁻¹, :⋆₁, :dual_d₁]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    Σ = 0
    sum = Int64[]
    Summand = 0
    summand = Int64[]
    summation = Int64[]
  end
  @test codiff_1 == test_codiff_1

  # Test for Codifferential (Primal2 -> Primal1)
  codiff_2 = @decapode begin
    A::Form2
    B::Form1

    B == δ₂(A)
  end
  Decapodes.open_operators!(codiff_2)
  infer_types!(codiff_2)
  resolve_overloads!(codiff_2)

  test_codiff_2 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 4
    name = [:A, :B, Symbol("•_1_1"), Symbol("•_1_2")]
    type = [:Form2, :Form1, :DualForm0, :DualForm1]
    Op1 = 3
    src = [4, 1, 3]
    tgt = [2, 3, 4]
    op1 = [:⋆₁⁻¹, :⋆₂, :dual_d₀]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    Σ = 0
    sum = Int64[]
    Summand = 0
    summand = Int64[]
    summation = Int64[]
  end
  @test codiff_2 == test_codiff_2

  # Test for Laplace de Rham (Primal0 -> Primal0)
  laplace_de_rham_0 = @decapode begin
    (A, B)::Form0

    B == Δ₀(A)
  end
  Decapodes.open_operators!(laplace_de_rham_0)
  infer_types!(laplace_de_rham_0)
  resolve_overloads!(laplace_de_rham_0)

  test_laplace_de_rham_0 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 5
    name = [:A, Symbol("•_1_1"), :B, Symbol("•_1_1"), Symbol("•_1_2")]
    type = [:Form0, :Form1, :Form0, :DualForm1, :DualForm2]
    Op1 = 4
    src = [5, 1, 2, 4]
    tgt = [3, 2, 4, 5]
    op1 = [:⋆₀⁻¹, :d₀, :⋆₁, :dual_d₁]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    Σ = 0
    sum = Int64[]
    Summand = 0
    summand = Int64[]
    summation = Int64[]
  end
  @test laplace_de_rham_0 == test_laplace_de_rham_0

  # Test for Laplace de Rham (Primal1 -> Primal1)
  laplace_de_rham_1 = @decapode begin
    (A, B)::Form1

    B == Δ₁(A)
  end
  Decapodes.open_operators!(laplace_de_rham_1)
  infer_types!(laplace_de_rham_1)
  resolve_overloads!(laplace_de_rham_1)

  test_laplace_de_rham_1 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 10
    name = [Symbol("•_2_1"), Symbol("•2"), :A, :B, Symbol("•_1_1"), Symbol("•1"), Symbol("•_2_1"), Symbol("•_2_2"), Symbol("•_3_1"), Symbol("•_3_2")]
    type = [:Form2, :Form1, :Form1, :Form1, :Form0, :Form1, :DualForm0, :DualForm1, :DualForm1, :DualForm2]
    Op1 = 8
    src = [3, 9, 10, 5, 1, 7, 8, 3]
    tgt = [1, 10, 5, 2, 7, 8, 6, 9]
    op1 = [:d₁, :dual_d₁, :⋆₀⁻¹, :d₀, :⋆₂, :dual_d₀, :⋆₁⁻¹, :⋆₁]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    Σ = 1
    sum = [4]
    Summand = 2
    summand = [6, 2]
    summation = [1, 1]
  end
  @test laplace_de_rham_1 == test_laplace_de_rham_1

  # Test for Laplace de Rham (Primal2 -> Primal2)
  laplace_de_rham_2 = @decapode begin
    (A, B)::Form2

    B == Δ₂(A)
  end
  Decapodes.open_operators!(laplace_de_rham_2)
  infer_types!(laplace_de_rham_2)
  resolve_overloads!(laplace_de_rham_2)

  test_laplace_de_rham_2 = @acset SummationDecapode{Any, Any, Symbol} begin
    Var = 5
    name = [:B, :A, Symbol("•_1_1"), Symbol("•_2_1"), Symbol("•_2_2")]
    type = [:Form2, :Form2, :Form1, :DualForm0, :DualForm1]    
    Op1 = 4
    src = [3, 5, 2, 4]
    tgt = [1, 3, 4, 5]
    op1 = [:d₁, :⋆₁⁻¹, :⋆₂, :dual_d₀]
    Op2 = 0
    proj1 = Int64[]
    proj2 = Int64[]
    res = Int64[]
    op2 = Any[]
    Σ = 0
    sum = Int64[]
    Summand = 0
    summand = Int64[]
    summation = Int64[]
  end
  @test laplace_de_rham_2 == test_laplace_de_rham_2
end
