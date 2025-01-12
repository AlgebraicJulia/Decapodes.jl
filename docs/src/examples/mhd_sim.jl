begin
    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:664 =#
    (mesh, operators, hodge = GeometricHodge())->begin
            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:664 =#
            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:665 =#
            begin
                #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:214 =#
                (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_matrix_generate(mesh, :⋆₀⁻¹, hodge)
                (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_matrix_generate(mesh, :dual_d₁, hodge)
                (var"GenSim-M_∧₀₁", ∧₀₁) = default_dec_matrix_generate(mesh, :∧₀₁, hodge)
                (∧ᵈᵈ₁₀) = operators(mesh, :∧ᵈᵈ₁₀)
                (var"GenSim-M_dual_d₀", dual_d₀) = default_dec_matrix_generate(mesh, :dual_d₀, hodge)
                (∧ᵖᵈ₁₁) = operators(mesh, :∧ᵖᵈ₁₁)
                (var"GenSim-M_⋆₁⁻¹", ⋆₁⁻¹) = default_dec_matrix_generate(mesh, :⋆₁⁻¹, hodge)
                (var"GenSim-M_⋆₁", ⋆₁) = default_dec_matrix_generate(mesh, :⋆₁, hodge)
                (var"GenSim-M_d₁", d₁) = default_dec_matrix_generate(mesh, :d₁, hodge)
                Δ⁻¹ = operators(mesh, :Δ⁻¹)
                (var"GenSim-M_⋆₂", ⋆₂) = default_dec_matrix_generate(mesh, :⋆₂, hodge)
                ♭♯ = operators(mesh, :♭♯)
                (var"GenSim-M_d₀", d₀) = default_dec_matrix_generate(mesh, :d₀, hodge)
            end
            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:666 =#
            begin
                #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:545 =#
                var"GenSim-M_GenSim-ConMat_0" = var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁"
                var"GenSim-ConMat_0" = (x->begin
                            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:562 =#
                            var"GenSim-M_GenSim-ConMat_0" * x
                        end)
                var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_dual_d₀" * var"GenSim-M_⋆₂"
                var"GenSim-ConMat_1" = (x->begin
                            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:562 =#
                            var"GenSim-M_GenSim-ConMat_1" * x
                        end)
            end
            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:667 =#
            begin
                #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:660 =#
            end
            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:668 =#
            begin
                #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:658 =#
                var"__•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__•_6_1" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__•_6_2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
                var"__•10" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
                var"__•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                var"__•_11_1" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
                var"__•12" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
                __η = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            end
            #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:669 =#
            f(__du__, __u__, __p__, __t__) = begin
                    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:669 =#
                    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:670 =#
                    begin
                        #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:253 =#
                        β = __u__.β
                        dη = __u__.dη
                        var"-1" = -1.0
                    end
                    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:671 =#
                    begin
                        #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:491 =#
                        var"•5" = Decapodes.get_tmp(var"__•5", __u__)
                        var"•9" = Decapodes.get_tmp(var"__•9", __u__)
                        var"•_6_1" = Decapodes.get_tmp(var"__•_6_1", __u__)
                        var"•_6_2" = Decapodes.get_tmp(var"__•_6_2", __u__)
                        var"•10" = Decapodes.get_tmp(var"__•10", __u__)
                        var"•2" = Decapodes.get_tmp(var"__•2", __u__)
                        var"•_11_1" = Decapodes.get_tmp(var"__•_11_1", __u__)
                        var"•12" = Decapodes.get_tmp(var"__•12", __u__)
                        η = Decapodes.get_tmp(__η, __u__)
                    end
                    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:672 =#
                    mul!(var"•5", var"GenSim-M_⋆₀⁻¹", dη)
                    var"GenSim-M_⋆₁⁻¹"(var"•9", β)
                    var"•8" = ♭♯(var"•9")
                    var"GenSim-M_⋆₁⁻¹"(var"•_6_1", β)
                    mul!(var"•_6_2", var"GenSim-M_d₁", var"•_6_1")
                    mul!(var"•10", var"GenSim-M_⋆₂", var"•_6_2")
                    var"GenSim-M_⋆₁⁻¹"(var"•2", β)
                    mul!(var"•_11_1", var"GenSim-M_⋆₀⁻¹", dη)
                    ψ = Δ⁻¹(var"•_11_1")
                    var"•7" = var"•8" ∧ᵈᵈ₁₀ var"•10"
                    mul!(var"•12", var"GenSim-M_d₀", ψ)
                    mul!(η, var"GenSim-M_⋆₁", var"•12")
                    var"•14" = var"•2" ∧ᵖᵈ₁₁ η
                    var"•6" = ♭♯(η)
                    var"•13" = var"GenSim-ConMat_1"(var"•14")
                    var"•4" = var"•5" ∧₀₁ var"•6"
                    β̇ = var"-1" .* var"•13"
                    sum_1 = (.+)(var"•4", var"•7")
                    var"•3" = var"GenSim-ConMat_0"(sum_1)
                    dη̇ = var"-1" .* var"•3"
                    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:673 =#
                    begin
                        #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:294 =#
                        setproperty!(__du__, :dη, dη̇)
                        setproperty!(__du__, :β, β̇)
                    end
                    #= /home/cuffaro/Documents/work/ufaj/dev/decapodes/docs-refactor/src/simulation.jl:674 =#
                    return nothing
                end
        end
end