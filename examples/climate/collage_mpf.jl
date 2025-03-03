begin
    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:692 =#
    (mesh, operators, hodge = GeometricHodge())->begin
            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:692 =#
            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:693 =#
            begin
                #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:217 =#
                (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_matrix_generate(mesh, :⋆₀⁻¹, hodge)
                Δ₀⁻¹ = default_dec_generate(mesh, :Δ₀⁻¹, hodge)
                (var"GenSim-M_d₀", d₀) = default_dec_matrix_generate(mesh, :d₀, hodge)
                (var"GenSim-M_⋆₁", ⋆₁) = default_dec_matrix_generate(mesh, :⋆₁, hodge)
                (var"GenSim-M_♭♯", ♭♯) = default_dec_matrix_generate(mesh, :♭♯, hodge)
                exp = operators(mesh, :exp)
                (∧ᵈᵖ₁₀) = default_dec_generate(mesh, :∧ᵈᵖ₁₀, hodge)
                bound_dual1form = operators(mesh, :bound_dual1form)
                bound_dual2form = operators(mesh, :bound_dual2form)
                (var"GenSim-M_∧₀₁", ∧₀₁) = default_dec_matrix_generate(mesh, :∧₀₁, hodge)
                (var"GenSim-M_d̃₁", d̃₁) = default_dec_matrix_generate(mesh, :d̃₁, hodge)
                (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_matrix_generate(mesh, :dual_d₁, hodge)
            end
            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:694 =#
            begin
                #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:552 =#
                var"GenSim-M_GenSim-ConMat_0" = var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁" * var"GenSim-M_d₀" * var"GenSim-M_⋆₀⁻¹"
                var"GenSim-ConMat_0" = (x->begin
                            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:569 =#
                            var"GenSim-M_GenSim-ConMat_0" * x
                        end)
                var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_d̃₁" * var"GenSim-M_⋆₁" * var"GenSim-M_♭♯"
                var"GenSim-ConMat_1" = (x->begin
                            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:569 =#
                            var"GenSim-M_GenSim-ConMat_1" * x
                        end)
                var"GenSim-M_GenSim-ConMat_2" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁" * var"GenSim-M_d₀"
                var"GenSim-ConMat_2" = (x->begin
                            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:569 =#
                            var"GenSim-M_GenSim-ConMat_2" * x
                        end)
                var"GenSim-M_GenSim-ConMat_3" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁"
                var"GenSim-ConMat_3" = (x->begin
                            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:569 =#
                            var"GenSim-M_GenSim-ConMat_3" * x
                        end)
            end
            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:695 =#
            begin
                #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:688 =#
            end
            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:696 =#
            begin
                #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:686 =#
                var"__phasefield_•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__phasefield_•6" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__phasefield_•8" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__phasefield_•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__navierstokes_•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__navierstokes_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                __𝐮 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__navierstokes_•8" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__phasefield_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__phasefield_•12" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__navierstokes_•6" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__phasefield_•3" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__phasefield_•11" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__phasefield_•10" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                __viscosity_sum_1 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                __phasefield_sum_1 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__phasefield_•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__viscosity_•1" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                __phasefield_Ċ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                __μ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
            end
            #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:697 =#
            f(__du__, __u__, __p__, __t__) = begin
                    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:697 =#
                    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:698 =#
                    begin
                        #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:261 =#
                        navierstokes_d𝐮 = __u__.navierstokes_d𝐮
                        navierstokes_U = __u__.navierstokes_U
                        navierstokes_DU = __u__.navierstokes_DU
                        C = __u__.C
                        viscosity_L = __p__.viscosity_L
                        viscosity_k = __p__.viscosity_k
                        viscosity_J = __p__.viscosity_J
                        phasefield_D = __p__.phasefield_D
                        phasefield_γ = __p__.phasefield_γ
                        phasefield_η = __p__.phasefield_η
                        phasefield_F = __p__.phasefield_F
                        var"1" = 1.0
                        var"3" = 3.0
                    end
                    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:699 =#
                    begin
                        #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:498 =#
                        var"phasefield_•9" = Decapodes.get_tmp(var"__phasefield_•9", __u__)
                        var"phasefield_•6" = Decapodes.get_tmp(var"__phasefield_•6", __u__)
                        var"phasefield_•8" = Decapodes.get_tmp(var"__phasefield_•8", __u__)
                        var"phasefield_•5" = Decapodes.get_tmp(var"__phasefield_•5", __u__)
                        var"navierstokes_•2" = Decapodes.get_tmp(var"__navierstokes_•2", __u__)
                        var"navierstokes_•4" = Decapodes.get_tmp(var"__navierstokes_•4", __u__)
                        𝐮 = Decapodes.get_tmp(__𝐮, __u__)
                        var"navierstokes_•8" = Decapodes.get_tmp(var"__navierstokes_•8", __u__)
                        var"phasefield_•4" = Decapodes.get_tmp(var"__phasefield_•4", __u__)
                        var"phasefield_•12" = Decapodes.get_tmp(var"__phasefield_•12", __u__)
                        var"navierstokes_•6" = Decapodes.get_tmp(var"__navierstokes_•6", __u__)
                        var"phasefield_•3" = Decapodes.get_tmp(var"__phasefield_•3", __u__)
                        var"phasefield_•11" = Decapodes.get_tmp(var"__phasefield_•11", __u__)
                        var"phasefield_•10" = Decapodes.get_tmp(var"__phasefield_•10", __u__)
                        viscosity_sum_1 = Decapodes.get_tmp(__viscosity_sum_1, __u__)
                        phasefield_sum_1 = Decapodes.get_tmp(__phasefield_sum_1, __u__)
                        var"phasefield_•2" = Decapodes.get_tmp(var"__phasefield_•2", __u__)
                        var"viscosity_•1" = Decapodes.get_tmp(var"__viscosity_•1", __u__)
                        phasefield_Ċ = Decapodes.get_tmp(__phasefield_Ċ, __u__)
                        μ = Decapodes.get_tmp(__μ, __u__)
                    end
                    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:700 =#
                    mul!(var"phasefield_•9", var"GenSim-M_GenSim-ConMat_2", C)
                    var"viscosity_•4" = -viscosity_k
                    navierstokes_r2_d𝐮 = bound_dual2form(navierstokes_d𝐮, navierstokes_DU)
                    var"viscosity_•3" = var"viscosity_•4" .* C
                    var"phasefield_•7" = C .^ var"3"
                    var"phasefield_•6" .= var"phasefield_•7" .- C
                    var"phasefield_•8" .= phasefield_γ .* var"phasefield_•9"
                    var"phasefield_•5" .= var"phasefield_•6" .- var"phasefield_•8"
                    mul!(var"navierstokes_•2", var"GenSim-M_⋆₀⁻¹", navierstokes_r2_d𝐮)
                    navierstokes_ψ = Δ₀⁻¹(var"navierstokes_•2")
                    mul!(var"navierstokes_•4", var"GenSim-M_d₀", navierstokes_ψ)
                    mul!(𝐮, var"GenSim-M_⋆₁", var"navierstokes_•4")
                    mul!(var"navierstokes_•8", var"GenSim-M_GenSim-ConMat_0", navierstokes_r2_d𝐮)
                    mul!(var"phasefield_•4", var"GenSim-M_d₀", var"phasefield_•5")
                    mul!(var"phasefield_•12", var"GenSim-M_♭♯", 𝐮)
                    mul!(var"navierstokes_•6", var"GenSim-M_⋆₀⁻¹", navierstokes_r2_d𝐮)
                    var"viscosity_•2" = exp(var"viscosity_•3")
                    navierstokes_r1_𝐮 = bound_dual1form(𝐮, navierstokes_U)
                    var"phasefield_•3" .= phasefield_F .* var"phasefield_•4"
                    var"GenSim-M_∧₀₁"(var"phasefield_•11", C, var"phasefield_•12")
                    var"phasefield_•10" .= phasefield_η .* var"phasefield_•11"
                    viscosity_sum_1 .= (.+)(var"1", var"viscosity_•2")
                    phasefield_sum_1 .= (.+)(var"phasefield_•3", var"phasefield_•10")
                    mul!(var"phasefield_•2", var"GenSim-M_GenSim-ConMat_3", phasefield_sum_1)
                    var"navierstokes_•3" = navierstokes_r1_𝐮 ∧ᵈᵖ₁₀ var"navierstokes_•6"
                    var"viscosity_•1" .= viscosity_L ./ viscosity_sum_1
                    phasefield_Ċ .= phasefield_D .* var"phasefield_•2"
                    μ .= (.+)(var"viscosity_•1", viscosity_J)
                    var"navierstokes_•1" = var"GenSim-ConMat_1"(var"navierstokes_•3")
                    var"navierstokes_•7" = μ .* var"navierstokes_•8"
                    navierstokes_d𝐮̇ = var"navierstokes_•7" .- var"navierstokes_•1"
                    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:701 =#
                    begin
                        #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:302 =#
                        setproperty!(__du__, :navierstokes_d𝐮, navierstokes_d𝐮̇)
                        setproperty!(__du__, :C, phasefield_Ċ)
                    end
                    #= /Users/hacker/Code/Decapodes.jl/src/simulation.jl:702 =#
                    return nothing
                end
        end
end