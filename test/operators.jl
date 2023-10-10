using Decapodes
using Test
using SparseArrays
using LinearAlgebra
using CombinatorialSpaces
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

import Decapodes: dec_p_differential, dec_p_laplace_de_rham, dec_p_hodge_diag, open_operators!

include("../examples/grid_meshes.jl")
include("../examples/sw/spherical_meshes.jl")

function generate_dual_mesh(s::HasDeltaSet2D)
    orient!(s)
    sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
    subdivide_duals!(sd, Barycenter())
    sd
end

dual_meshes = [(generate_dual_mesh ∘ loadmesh ∘ Icosphere).(1:3)...,
               (generate_dual_mesh ∘ loadmesh)(Rectangle_30x10()),
               (generate_dual_mesh).([triangulated_grid(10,10,8,8,Point3D), makeSphere(5, 175, 5, 0, 360, 5, 6371+90)[1]])...,
               (loadmesh)(Torus_30x10())]

@testset "In-House Exterior Derivative for Form0" begin    
    for sd in dual_meshes
        @test dec_p_differential(Val{0}, sd) == d(0, sd)
    end
end

@testset "In-House Exterior Derivative for Form1" begin
    for sd in dual_meshes
        @test dec_p_differential(Val{1}, sd) == d(1, sd)
    end
end

@testset "In-House Laplace DeRham for Form1" begin
    for sd in dual_meshes
        @test dec_p_laplace_de_rham(0, sd) == Δ(0, sd)
    end
end

#TODO: The laplace DeRham on Torus produces a non-singular matrix
@testset "In-House Laplace DeRham for Form1" begin
    for sd in dual_meshes[1:6]
        @test dec_p_laplace_de_rham(1, sd) == Δ(1, sd)
    end
end

@testset "In-House Diagonal Hodge for Form1" begin
    for sd in dual_meshes
        @test dec_p_hodge_diag(Val{1}, sd) == hodge_star(1, sd, DiagonalHodge())
    end
end

@testset "In-House Diagonal Hodge for Form2" begin
    for sd in dual_meshes
        @test dec_p_hodge_diag(Val{2}, sd) == hodge_star(2, sd, DiagonalHodge())
    end
end

function print_decapode_into_acset(d::SummationDecapode)
    println("@acset SummationDecapode{Any, Any, Symbol} begin")
    println("   Var = ", nparts(d, :Var))
    println("   name = ", d[:name])
    println("   type = ", d[:type])

    println("   Op1 = ", nparts(d, :Op1))
    println("   src = ", d[:src])
    println("   tgt = ", d[:tgt])
    println("   op1 = ", d[:op1])

    println("   Op2 = ", nparts(d, :Op2))
    println("   proj1 = ", d[:proj1])
    println("   proj2 = ", d[:proj2])
    println("   res = ", d[:res])
    println("   op2 = ", d[:op2])

    println("   Σ = ", nparts(d, :Σ))
    println("   sum = ", d[:sum])

    println("   Summand = ", nparts(d, :Summand))
    println("   summand = ", d[:summand])
    println("   summation = ", d[:summation])
    println("end")
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

    lie_derivative_1 = @decapode begin
        A::DualForm1
        B::Form1
        C::DualForm1

        C == L₁(B, A)
    end
    open_operators!(lie_derivative_1)
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

    lie_derivative_2 = @decapode begin
        A::DualForm2
        B::Form1
        C::DualForm2

        C == L₂(B, A)
    end
    open_operators!(lie_derivative_2)
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

end