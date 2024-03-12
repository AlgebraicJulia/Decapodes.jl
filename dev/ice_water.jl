begin
    #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:531 =#
    function simulate(mesh, operators, hodge = GeometricHodge())
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:531 =#
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:532 =#
        begin
            #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:157 =#
            (var"GenSim-M_d₀", d₀) = default_dec_matrix_generate(mesh, :d₀, hodge)
            (var"GenSim-M_⋆₁", ⋆₁) = default_dec_matrix_generate(mesh, :⋆₁, hodge)
            (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_matrix_generate(mesh, :dual_d₁, hodge)
            (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_matrix_generate(mesh, :⋆₀⁻¹, hodge)
            (var"GenSim-M_⋆₁⁻¹", ⋆₁⁻¹) = default_dec_matrix_generate(mesh, :⋆₁⁻¹, hodge)
            (var"GenSim-M_dual_d₀", dual_d₀) = default_dec_matrix_generate(mesh, :dual_d₀, hodge)
            (var"GenSim-M_∧₁₀", ∧₁₀) = default_dec_matrix_generate(mesh, :∧₁₀, hodge)
            ♯ = operators(mesh, :♯)
            mag = operators(mesh, :mag)
            σ = operators(mesh, :σ)
            (^) = operators(mesh, :^)
            ι₁₁ = operators(mesh, :ι₁₁)
            (∧ᵈᵖ₁₀) = operators(mesh, :∧ᵈᵖ₁₀)
            (∧ᵖᵈ₀₁) = operators(mesh, :∧ᵖᵈ₀₁)
        end
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:533 =#
        begin
            #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:449 =#
            var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁"
            var"GenSim-ConMat_1" = (x->var"GenSim-M_GenSim-ConMat_1" * x)
            var"GenSim-M_GenSim-ConMat_2" = var"GenSim-M_⋆₁" * var"GenSim-M_d₀" * var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁"
            var"GenSim-ConMat_2" = (x->var"GenSim-M_GenSim-ConMat_2" * x)
            var"GenSim-M_GenSim-ConMat_3" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁"
            var"GenSim-ConMat_3" = (x->var"GenSim-M_GenSim-ConMat_3" * x)
        end
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:534 =#
        begin
            #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:203 =#
            var"__glacier_dynamics_dynamics_•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__glacier_dynamics_dynamics_•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__water_dynamics_•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
            var"__water_dynamics_•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__water_dynamics_•1" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__water_dynamics_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __water_dynamics_𝑝ᵈ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__water_dynamics_•7" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __water_dynamics_sum_1 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__glacier_dynamics_dynamics_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__glacier_dynamics_dynamics_•3" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__glacier_dynamics_dynamics_•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__water_dynamics_•6" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __water_dynamics_𝐮̇ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __glacier_dynamics_dynamics_ḣ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
        end
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:535 =#
        f(du, u, p, t) = begin
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:535 =#
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:536 =#
                begin
                    #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:228 =#
                    ice_thickness = u.ice_thickness
                    glacier_dynamics_n = p.glacier_dynamics_n
                    glacier_dynamics_stress_A = p.glacier_dynamics_stress_A
                    glacier_dynamics_stress_ρ = p.glacier_dynamics_stress_ρ
                    glacier_dynamics_stress_g = p.glacier_dynamics_stress_g
                    flow = u.flow
                    water_dynamics_P = u.water_dynamics_P
                    water_dynamics_μ = p.water_dynamics_μ
                    var"1" = 1.0
                    var"2" = 2.0
                    var"2" = 2.0
                    var"0.5" = 0.5
                    var"-1" = -1.0
                    var"1" = 1.0
                end
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:537 =#
                var"glacier_dynamics_dynamics_•5" = Decapodes.get_tmp(var"__glacier_dynamics_dynamics_•5", u)
                var"glacier_dynamics_dynamics_•9" = Decapodes.get_tmp(var"__glacier_dynamics_dynamics_•9", u)
                var"water_dynamics_•9" = Decapodes.get_tmp(var"__water_dynamics_•9", u)
                var"water_dynamics_•5" = Decapodes.get_tmp(var"__water_dynamics_•5", u)
                var"water_dynamics_•1" = Decapodes.get_tmp(var"__water_dynamics_•1", u)
                var"water_dynamics_•4" = Decapodes.get_tmp(var"__water_dynamics_•4", u)
                water_dynamics_𝑝ᵈ = Decapodes.get_tmp(__water_dynamics_𝑝ᵈ, u)
                var"water_dynamics_•7" = Decapodes.get_tmp(var"__water_dynamics_•7", u)
                water_dynamics_sum_1 = Decapodes.get_tmp(__water_dynamics_sum_1, u)
                var"glacier_dynamics_dynamics_•4" = Decapodes.get_tmp(var"__glacier_dynamics_dynamics_•4", u)
                var"glacier_dynamics_dynamics_•3" = Decapodes.get_tmp(var"__glacier_dynamics_dynamics_•3", u)
                var"glacier_dynamics_dynamics_•2" = Decapodes.get_tmp(var"__glacier_dynamics_dynamics_•2", u)
                var"water_dynamics_•6" = Decapodes.get_tmp(var"__water_dynamics_•6", u)
                water_dynamics_𝐮̇ = Decapodes.get_tmp(__water_dynamics_𝐮̇, u)
                glacier_dynamics_dynamics_ḣ = Decapodes.get_tmp(__glacier_dynamics_dynamics_ḣ, u)
                mul!(var"glacier_dynamics_dynamics_•5", var"GenSim-M_d₀", ice_thickness)
                mul!(var"glacier_dynamics_dynamics_•9", var"GenSim-M_d₀", ice_thickness)
                var"glacier_dynamics_dynamics_•8" = ♯(var"glacier_dynamics_dynamics_•9")
                var"glacier_dynamics_dynamics_•7" = mag(var"glacier_dynamics_dynamics_•8")
                var"interaction_•1" = σ(ice_thickness)
                var"glacier_dynamics_dynamics_•10" = glacier_dynamics_n .- var"1"
                var"glacier_dynamics_dynamics_•6" = var"glacier_dynamics_dynamics_•7" ^ var"glacier_dynamics_dynamics_•10"
                var"glacier_dynamics_stress_•3" = glacier_dynamics_stress_ρ .* glacier_dynamics_stress_g
                var"glacier_dynamics_stress_•2" = var"glacier_dynamics_stress_•3" ^ glacier_dynamics_n
                var"interaction_•2" = var"1" .- var"interaction_•1"
                flow_after = var"interaction_•2" ∧ᵖᵈ₀₁ flow
                glacier_dynamics_dynamics_sum_1 = (.+)(glacier_dynamics_n, var"2")
                glacier_dynamics_stress_sum_1 = (.+)(glacier_dynamics_n, var"2")
                mul!(var"water_dynamics_•9", var"GenSim-M_GenSim-ConMat_1", flow_after)
                mul!(var"water_dynamics_•5", var"GenSim-M_GenSim-ConMat_2", flow_after)
                var"glacier_dynamics_dynamics_•11" = ice_thickness ^ glacier_dynamics_dynamics_sum_1
                var"glacier_dynamics_stress_•1" = var"2" / glacier_dynamics_stress_sum_1
                glacier_dynamics_stress_mult_1 = var"glacier_dynamics_stress_•1" .* glacier_dynamics_stress_A
                glacier_dynamics_Γ = glacier_dynamics_stress_mult_1 .* var"glacier_dynamics_stress_•2"
                var"water_dynamics_•2" = ι₁₁(flow_after, flow_after)
                var"water_dynamics_•1" .= var"0.5" .* var"water_dynamics_•2"
                var"water_dynamics_•4" .= water_dynamics_μ .* var"water_dynamics_•5"
                var"water_dynamics_•8" = flow_after ∧ᵈᵖ₁₀ var"water_dynamics_•9"
                water_dynamics_𝑝ᵈ .= (.+)(water_dynamics_P, var"water_dynamics_•1")
                var"GenSim-M_⋆₁⁻¹"(var"water_dynamics_•7", var"water_dynamics_•8")
                mul!(water_dynamics_sum_1, var"GenSim-M_dual_d₀", water_dynamics_𝑝ᵈ)
                var"glacier_dynamics_dynamics_•4" .= glacier_dynamics_Γ .* var"glacier_dynamics_dynamics_•5"
                var"GenSim-M_∧₁₀"(var"glacier_dynamics_dynamics_•3", var"glacier_dynamics_dynamics_•4", var"glacier_dynamics_dynamics_•6")
                var"GenSim-M_∧₁₀"(var"glacier_dynamics_dynamics_•2", var"glacier_dynamics_dynamics_•3", var"glacier_dynamics_dynamics_•11")
                var"water_dynamics_•6" .= var"-1" .* var"water_dynamics_•7"
                water_dynamics_𝐮̇ .= (.+)(var"water_dynamics_•4", var"water_dynamics_•6", water_dynamics_sum_1)
                mul!(glacier_dynamics_dynamics_ḣ, var"GenSim-M_GenSim-ConMat_3", var"glacier_dynamics_dynamics_•2")
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:538 =#
                getproperty(du, :ice_thickness) .= glacier_dynamics_dynamics_ḣ
                getproperty(du, :flow) .= water_dynamics_𝐮̇
            end
    end
end