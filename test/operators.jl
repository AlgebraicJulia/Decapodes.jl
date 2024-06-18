@testset "Open Operators" begin
    laplace_de_rham_0 = @decapode begin
        (A, B)::Form0

        B == Δ₀(A)
    end
    test_compare = open_operators(laplace_de_rham_0)
    infer_types!(test_compare)
    resolve_overloads!(test_compare)

    open_operators!(laplace_de_rham_0)
    infer_types!(laplace_de_rham_0)
    resolve_overloads!(laplace_de_rham_0)


    @test test_compare == laplace_de_rham_0
    @test test_compare !== laplace_de_rham_0

end

@testset "Opening 1D Operators" begin
    # Test for 1D Interior Product (Dual1, Primal1 -> Dual0)
    interior_product_1 = @decapode begin
        A::DualForm1
        B::Form1
        C::DualForm0

        C == i₁(A, B)
    end
    open_operators!(interior_product_1; dimension = 1)
    infer_types!(interior_product_1, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(interior_product_1, op1_res_rules_1D, op2_res_rules_1D)

    test_interior_product_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 5
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2]
        type = [:DualForm1, :Form1, :DualForm0, :Form0, :Form1]
        Op1 = 2
        src = [1, 5]
        tgt = [4, 3]
        op1 = [:⋆₀⁻¹, :⋆₁]
        Op2 = 1
        proj1 = [4]
        proj2 = [2]
        res = [5]
        op2 = [:∧₀₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test interior_product_1 == test_interior_product_1

    # Test for 1D Lie Derivative (Dual0, Primal1 -> Dual0)
    lie_derivative_0 = @decapode begin
        A::DualForm0
        B::Form1
        C::DualForm0

        C == L₀(B, A)
    end
    open_operators!(lie_derivative_0; dimension = 1)
    infer_types!(lie_derivative_0, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(lie_derivative_0, op1_res_rules_1D, op2_res_rules_1D)

    test_lie_derivative_0 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 6
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]
        type = [:DualForm0, :Form1, :DualForm0, :DualForm1, :Form0, :Form1]
        Op1 = 3
        src = [1, 4, 6]
        tgt = [4, 5, 3]
        op1 = [:dual_d₀, :⋆₀⁻¹, :⋆₁]
        Op2 = 1
        proj1 = [5]
        proj2 = [2]
        res = [6]
        op2 = [:∧₀₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test lie_derivative_0 == test_lie_derivative_0

    # Test for 1D Lie Derivative (Dual1, Primal1 -> Dual1)
    lie_derivative_1 = @decapode begin
        A::DualForm1
        B::Form1
        C::DualForm1

        C == L₁(B, A)
    end
    open_operators!(lie_derivative_1; dimension = 1)
    infer_types!(lie_derivative_1, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(lie_derivative_1, op1_res_rules_1D, op2_res_rules_1D)

    test_lie_derivative_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 6
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]
        type = [:DualForm1, :Form1, :DualForm1, :DualForm0, :Form0, :Form1]
        Op1 = 3
        src = [1, 6, 4]
        tgt = [5, 4, 3]
        op1 = [:⋆₀⁻¹, :⋆₁, :dual_d₀]
        Op2 = 1
        proj1 = [5]
        proj2 = [2]
        res = [6]
        op2 = [:∧₀₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test lie_derivative_1 == test_lie_derivative_1

    # Test for Codifferential (Primal1 -> Primal0)
    codiff_1 = @decapode begin
        A::Form1
        B::Form0

        B == δ₁(A)
    end
    open_operators!(codiff_1, dimension = 1)
    infer_types!(codiff_1, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(codiff_1, op1_res_rules_1D, op2_res_rules_1D)

    test_codiff_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 4
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2]
        type = [:Form1, :Form0, :DualForm0, :DualForm1]
        Op1 = 3
        src = [4, 1, 3]
        tgt = [2, 3, 4]
        op1 = [:⋆₀⁻¹, :⋆₁, :dual_d₀]
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
    
    # Test for Laplace de Rham (Primal0 -> Primal0)
    laplace_de_rham_0 = @decapode begin
        (A, B)::Form0

        B == Δ₀(A)
    end
    open_operators!(laplace_de_rham_0, dimension = 1)
    infer_types!(laplace_de_rham_0, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(laplace_de_rham_0, op1_res_rules_1D, op2_res_rules_1D)

    test_laplace_de_rham_0 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 5
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]
        type = [:Form0, :Form0, :Form1, :DualForm0, :DualForm1]        
        Op1 = 4
        src = [5, 1, 3, 4]
        tgt = [2, 3, 4, 5]
        op1 = [:⋆₀⁻¹, :d₀, :⋆₁, :dual_d₀]
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
    open_operators!(laplace_de_rham_1, dimension = 1)
    infer_types!(laplace_de_rham_1, op1_inf_rules_1D, op2_inf_rules_1D)
    resolve_overloads!(laplace_de_rham_1, op1_res_rules_1D, op2_res_rules_1D)

    test_laplace_de_rham_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 5
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]   
        type = [:Form1, :Form1, :Form0, :DualForm0, :DualForm1]        
        Op1 = 4
        src = [3, 1, 4, 5]
        tgt = [2, 4, 5, 3]
        op1 = [:d₀, :⋆₁, :dual_d₀, :⋆₀⁻¹]
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
    @test laplace_de_rham_1 == test_laplace_de_rham_1
end

@testset "Opening 2D Operators" begin
    # Test for Interior Product (Dual1, Primal1 -> Dual0)
    interior_product_1 = @decapode begin
        A::DualForm1
        B::Form1
        C::DualForm0

        C == i₁(A, B)
    end
    open_operators!(interior_product_1)
    infer_types!(interior_product_1)
    resolve_overloads!(interior_product_1)

    test_interior_product_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 6
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]
        type = [:DualForm1, :Form1, :DualForm0, :DualForm0, :Form1, :Form2]
        Op1 = 3
        src = [1, 6, 4]
        tgt = [5, 4, 3]
        op1 = [:⋆₁⁻¹, :⋆₂, :neg]
        Op2 = 1
        proj1 = [5]
        proj2 = [2]
        res = [6]
        op2 = [:∧₁₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test interior_product_1 == test_interior_product_1

    # Test for Interior Product (Dual2, Primal1 -> Dual1)
    interior_product_2 = @decapode begin
        A::DualForm2
        B::Form1
        C::DualForm1

        C == i₂(A, B)
    end
    open_operators!(interior_product_2)
    infer_types!(interior_product_2)
    resolve_overloads!(interior_product_2)

    test_interior_product_2 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 5
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2]
        type = [:DualForm2, :Form1, :DualForm1, :Form0, :Form1]
        Op1 = 2
        src = [1, 5]
        tgt = [4, 3]
        op1 = [:⋆₀⁻¹, :⋆₁]
        Op2 = 1
        proj1 = [4]
        proj2 = [2]
        res = [5]
        op2 = [:∧₀₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test interior_product_2 == test_interior_product_2

    # Test for Lie Derivative (Primal1, Dual0 -> Dual0)
    lie_derivative_0 = @decapode begin
        A::DualForm0
        B::Form1
        C::DualForm0

        C == L₀(B, A)
    end
    open_operators!(lie_derivative_0)
    infer_types!(lie_derivative_0)
    resolve_overloads!(lie_derivative_0)

    test_lie_derivative_0 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 6
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]
        type = [:DualForm0, :Form1, :DualForm0, :DualForm1, :Form1, :Form2]
        Op1 = 3
        src = [1, 4, 6]
        tgt = [4, 5, 3]
        op1 = [:dual_d₀, :⋆₁⁻¹, :⋆₂]
        Op2 = 1
        proj1 = [5]
        proj2 = [2]
        res = [6]
        op2 = [:∧₁₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test lie_derivative_0 == test_lie_derivative_0

    # Test for Lie Derivative (Primal1, Dual1 -> Dual1)
    lie_derivative_1 = @decapode begin
        A::DualForm1
        B::Form1
        C::DualForm1

        C == L₁(B, A)
    end
    open_operators!(lie_derivative_1)
    infer_types!(lie_derivative_1)
    resolve_overloads!(lie_derivative_1)

    test_lie_derivative_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 12
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3, :Gensim_Var_4, :Gensim_Var_5, :Gensim_Var_6, :Gensim_Var_7, :Gensim_Var_8, :Gensim_Var_9]
        type = [:DualForm1, :Form1, :DualForm1, :DualForm2, :DualForm1, :Form0, :Form1, :DualForm0, :DualForm0, :Form1, :Form2, :DualForm1]
        Op1 = 7
        src = [1, 4, 7, 1, 11, 9, 8]
        tgt = [4, 6, 5, 10, 9, 8, 12]
        op1 = [:dual_d₁, :⋆₀⁻¹, :⋆₁, :⋆₁⁻¹, :⋆₂, :neg, :dual_d₀]
        Op2 = 2
        proj1 = [10, 6]
        proj2 = [2, 2]
        res = [11, 7]
        op2 = [:∧₁₁, :∧₀₁]
        Σ = 1
        sum = [3]
        Summand = 2
        summand = [5, 12]
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
    open_operators!(lie_derivative_2)
    infer_types!(lie_derivative_2)
    resolve_overloads!(lie_derivative_2)

    test_lie_derivative_2 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 6
        name = [:A, :B, :C, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]
        type = [:DualForm2, :Form1, :DualForm2, :DualForm1, :Form0, :Form1]
        Op1 = 3
        src = [1, 6, 4]
        tgt = [5, 4, 3]
        op1 = [:⋆₀⁻¹, :⋆₁, :dual_d₁]
        Op2 = 1
        proj1 = [5]
        proj2 = [2]
        res = [6]
        op2 = [:∧₀₁]
        Σ = 0
        sum = Int64[]
        Summand = 0
        summand = Int64[]
        summation = Int64[]
    end
    @test lie_derivative_2 == test_lie_derivative_2

    # Test for Codifferential (Primal1 -> Primal0)
    codiff_1 = @decapode begin
        A::Form1
        B::Form0

        B == δ₁(A)
    end
    open_operators!(codiff_1)
    infer_types!(codiff_1)
    resolve_overloads!(codiff_1)

    test_codiff_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 4
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2]
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
    open_operators!(codiff_2)
    infer_types!(codiff_2)
    resolve_overloads!(codiff_2)

    test_codiff_2 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 4
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2]
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
    open_operators!(laplace_de_rham_0)
    infer_types!(laplace_de_rham_0)
    resolve_overloads!(laplace_de_rham_0)

    test_laplace_de_rham_0 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 5
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]   
        type = [:Form0, :Form0, :Form1, :DualForm1, :DualForm2]        
        Op1 = 4
        src = [5, 1, 3, 4]
        tgt = [2, 3, 4, 5]
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
    open_operators!(laplace_de_rham_1)
    infer_types!(laplace_de_rham_1)
    resolve_overloads!(laplace_de_rham_1)

    test_laplace_de_rham_1 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 10
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3, :Gensim_Var_4, :Gensim_Var_5, :Gensim_Var_6, :Gensim_Var_7, :Gensim_Var_8]
        type = [:Form1, :Form1, :Form1, :Form2, :DualForm0, :DualForm1, :Form1, :Form0, :DualForm1, :DualForm2]
        Op1 = 8
        src = [8, 1, 4, 5, 6, 1, 9, 10]
        tgt = [7, 4, 5, 6, 3, 9, 10, 8]
        op1 = [:d₀, :d₁, :⋆₂, :dual_d₀, :⋆₁⁻¹, :⋆₁, :dual_d₁, :⋆₀⁻¹]   
        Op2 = 0
        proj1 = Int64[]
        proj2 = Int64[]
        res = Int64[]
        op2 = Any[]
        Σ = 1
        sum = [2]
        Summand = 2
        summand = [3, 7]
        summation = [1, 1]
    end
    @test laplace_de_rham_1 == test_laplace_de_rham_1

    # Test for Laplace de Rham (Primal2 -> Primal2)
    laplace_de_rham_2 = @decapode begin
        (A, B)::Form2

        B == Δ₂(A)
    end
    open_operators!(laplace_de_rham_2)
    infer_types!(laplace_de_rham_2)
    resolve_overloads!(laplace_de_rham_2)

    test_laplace_de_rham_2 = @acset SummationDecapode{Any, Any, Symbol} begin
        Var = 5
        name = [:A, :B, :Gensim_Var_1, :Gensim_Var_2, :Gensim_Var_3]   
        type = [:Form2, :Form2, :Form1, :DualForm0, :DualForm1]        
        Op1 = 4
        src = [3, 1, 4, 5]
        tgt = [2, 4, 5, 3]
        op1 = [:d₁, :⋆₂, :dual_d₀, :⋆₁⁻¹]
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
