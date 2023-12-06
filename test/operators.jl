using Decapodes
using Test
using SparseArrays
using LinearAlgebra
using CombinatorialSpaces
using Random
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

import Decapodes: dec_wedge_product, dec_differential, dec_boundary, dec_dual_derivative, dec_hodge_star, dec_inv_hodge, open_operators!, open_operators
include("../examples/grid_meshes.jl")
include("../examples/sw/spherical_meshes.jl")

Random.seed!(0)

function generate_dual_mesh(s::HasDeltaSet1D)
    orient!(s)
    sd = EmbeddedDeltaDualComplex1D{Bool,Float64,Point2D}(s)
    subdivide_duals!(sd, Barycenter())
    sd
end

function generate_dual_mesh(s::HasDeltaSet2D)
    orient!(s)
    sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
    subdivide_duals!(sd, Barycenter())
    sd
end

primal_line = EmbeddedDeltaSet1D{Bool,Point2D}()
add_vertices!(primal_line, 3, point=[Point2D(1,0), Point2D(0,0), Point2D(0,2)])
add_edges!(primal_line, [1,2], [2,3])
line = generate_dual_mesh(primal_line)

primal_cycle = EmbeddedDeltaSet1D{Bool,Point2D}()
add_vertices!(primal_cycle, 3, point=[Point2D(1,0), Point2D(0,0), Point2D(0,1)])
add_edges!(primal_cycle, [1,2,3], [2,3,1])
cycle = generate_dual_mesh(primal_cycle)

primal_plus = EmbeddedDeltaSet1D{Bool,Point2D}()
add_vertices!(primal_plus, 5, point=[Point2D(0,0), Point2D(1,0), Point2D(-1,0), Point2D(0,1), Point2D(0, -1)])
add_edges!(primal_plus, [1,1,3,5], [2,4,1,1])
primal_plus[:edge_orientation] = true
plus = generate_dual_mesh(primal_plus)


dual_meshes_1D = [line, cycle, plus]

dual_meshes_2D = [(generate_dual_mesh ∘ loadmesh ∘ Icosphere).(1:2)...,
               (generate_dual_mesh ∘ loadmesh)(Rectangle_30x10()),
               (generate_dual_mesh).([triangulated_grid(10,10,8,8,Point3D), makeSphere(5, 175, 5, 0, 360, 5, 6371+90)[1]])...,
               (loadmesh)(Torus_30x10())];

@testset "Exterior Derivative" begin
    for i in 0:0 
        for sd in dual_meshes_1D
            @test all(dec_differential(i, sd) .== d(i, sd))
        end
    end

    for i in 0:1 
        for sd in dual_meshes_2D
            @test all(dec_differential(i, sd) .== d(i, sd))
        end
    end
end

@testset "Boundary" begin
    for i in 1:1 
        for sd in dual_meshes_1D
            @test all(dec_boundary(i, sd) .== ∂(i, sd))
        end
    end

    for i in 1:2
        for sd in dual_meshes_2D
            @test all(dec_boundary(i, sd) .== ∂(i, sd))
        end
    end
end

@testset "Dual Derivative" begin
    for i in 0:0 
        for sd in dual_meshes_1D
            @test all(dec_dual_derivative(i, sd) .== dual_derivative(i, sd))
        end
    end

    for i in 0:1
        for sd in dual_meshes_2D
            @test all(dec_dual_derivative(i, sd) .== dual_derivative(i, sd))
        end
    end
end

#TODO: For hodge star 1, the values seems to be extremely close yet not quite equal
@testset "Diagonal Hodge" begin
    for i in 0:1
        for sd in dual_meshes_1D
            @test all(isapprox.(dec_hodge_star(Val{i}, sd, DiagonalHodge()), hodge_star(i, sd, DiagonalHodge()); rtol = 1e-15))
        end
    end

    for i in 0:2
        for sd in dual_meshes_2D
            @test all(isapprox.(dec_hodge_star(Val{i}, sd, DiagonalHodge()), hodge_star(i, sd, DiagonalHodge()); rtol = 1e-15))
        end
    end
end

#TODO: For inv hodge star 1, the values seems to be extremely close yet not quite equal
@testset "Inverse Diagonal Hodge" begin
    for i in 0:1
        for sd in dual_meshes_1D
            @test all(isapprox.(dec_inv_hodge(Val{i}, sd, DiagonalHodge()), inv_hodge_star(i, sd, DiagonalHodge()); rtol = 1e-15))
        end
    end

    for i in 0:2
        for sd in dual_meshes_2D
            @test all(isapprox.(dec_inv_hodge(Val{i}, sd, DiagonalHodge()), inv_hodge_star(i, sd, DiagonalHodge()); rtol = 1e-15))
        end
    end
end

@testset "Geometric Hodge" begin
    for i in 0:1
        for sd in dual_meshes_1D
            @test all(isapprox.(dec_hodge_star(Val{i}, sd, GeometricHodge()), hodge_star(i, sd, GeometricHodge()); rtol = 1e-15))
        end
    end

    for i in [0, 2]
        for sd in dual_meshes_2D
            @test all(isapprox.(dec_hodge_star(Val{i}, sd, GeometricHodge()), hodge_star(i, sd, GeometricHodge()); rtol = 1e-15))
        end
    end

    # TODO: Why does this test require atol, not rtol, to reasonably pass?
    for i in [1]
        for sd in dual_meshes_2D
            @test all(isapprox.(dec_hodge_star(Val{i}, sd, GeometricHodge()), hodge_star(i, sd, GeometricHodge()); atol = 1e-15))
        end
    end

end

@testset "Inverse Geometric Hodge" begin
    for i in 0:1
        for sd in dual_meshes_1D
            @test all(isapprox.(dec_inv_hodge(Val{i}, sd, GeometricHodge()), inv_hodge_star(i, sd, GeometricHodge()); rtol = 1e-15))
        end
    end

    for i in [0, 2]
        for sd in dual_meshes_2D
            @test all(isapprox.(dec_inv_hodge(Val{i}, sd, GeometricHodge()), inv_hodge_star(i, sd, GeometricHodge()); rtol = 1e-15))
        end
    end

    for i in 1:1
        for sd in dual_meshes_2D[1:end-1]
            @show V_1 = rand(ne(sd))
            @test all(isapprox.(dec_inv_hodge(Val{i}, sd, GeometricHodge())(V_1), inv_hodge_star(i, sd, GeometricHodge()) * V_1; rtol = 1e-13))
        end
    end

end

@testset "Wedge Product" begin
    for sd in dual_meshes_1D
        V_1, V_2 = rand(nv(sd)), rand(nv(sd))
        E_1 = rand(ne(sd))
        @test all(dec_wedge_product(Tuple{0, 0}, sd)(V_1, V_2) .== ∧(Tuple{0, 0}, sd, V_1, V_2))
        @test all(dec_wedge_product(Tuple{0, 1}, sd)(V_1, E_1) .== ∧(Tuple{0, 1}, sd, V_1, E_1))
    end

    for sd in dual_meshes_2D
        V_1, V_2 = rand(nv(sd)), rand(nv(sd))
        E_1, E_2 = rand(ne(sd)), rand(ne(sd))
        T_2 = rand(ntriangles(sd))
        V_ones = ones(nv(sd))
        E_ones = ones(ne(sd))
        @test all(dec_wedge_product(Tuple{0, 0}, sd)(V_1, V_2) .== ∧(Tuple{0, 0}, sd, V_1, V_2))

        wdg01 = dec_wedge_product(Tuple{0, 1}, sd)
        @test all(isapprox.(wdg01(V_1, E_2), ∧(Tuple{0, 1}, sd, V_1, E_2); rtol = 1e-14))
        @test all(wdg01(V_ones, E_ones) .== E_ones)

        @test all(dec_wedge_product(Tuple{0, 2}, sd)(V_1, T_2) / 2 .== ∧(Tuple{0, 2}, sd, V_1, T_2))

        @test all(dec_wedge_product(Tuple{1, 1}, sd)(E_1, E_2) / 2 .== ∧(Tuple{1, 1}, sd, E_1, E_2))
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