begin
    #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:531 =#
    function simulate(mesh, operators, hodge = GeometricHodge())
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:531 =#
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:532 =#
        begin
            #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:157 =#
            (var"GenSim-M_dual_d₀", dual_d₀) = default_dec_matrix_generate(mesh, :dual_d₀, hodge)
            (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_matrix_generate(mesh, :dual_d₁, hodge)
            (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_matrix_generate(mesh, :⋆₀⁻¹, hodge)
            Δᵈ₁ = operators(mesh, :Δᵈ₁)
            Δᵈ₀ = operators(mesh, :Δᵈ₀)
            ℒ₁ = operators(mesh, :ℒ₁)
            ι₁₁ = operators(mesh, :ι₁₁)
            ι₁₂ = operators(mesh, :ι₁₂)
            (∧ᵖᵈ₀₁) = operators(mesh, :∧ᵖᵈ₀₁)
            (∧ᵈᵈ₀₁) = operators(mesh, :∧ᵈᵈ₀₁)
        end
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:533 =#
        begin
            #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:449 =#
            var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁"
            var"GenSim-ConMat_1" = (x->var"GenSim-M_GenSim-ConMat_1" * x)
        end
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:534 =#
        begin
            #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:203 =#
            var"__salinity_continuity_•11" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•14" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
            var"__momentum_•16" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
            var"__salinity_continuity_•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•19" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
            var"__momentum_•20" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__temperature_continuity_•7" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__temperature_continuity_•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__temperature_continuity_•11" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__salinity_continuity_•7" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•18" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :V)))
            __SD = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __temperature_FD = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__temperature_continuity_•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__temperature_continuity_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__temperature_continuity_•3" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__temperature_continuity_•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            __SD = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __salinity_FD = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__salinity_continuity_•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__salinity_continuity_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__salinity_continuity_•3" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__salinity_continuity_•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__momentum_•6" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __temperature_continuity_ċ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            __salinity_continuity_ċ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :Tri)))
            var"__momentum_•9" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•11" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•8" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __momentum_sum_1 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•5" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __momentum_sum_2 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•4" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•3" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __momentum_sum_3 = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            var"__momentum_•2" = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
            __momentum_v̇ = Decapodes.FixedSizeDiffCache(Vector{Float64}(undef, nparts(mesh, :E)))
        end
        #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:535 =#
        f(du, u, p, t) = begin
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:535 =#
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:536 =#
                begin
                    #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:228 =#
                    v = u.v
                    V = u.V
                    momentum_f = u.momentum_f
                    momentum_uˢ = u.momentum_uˢ
                    momentum_∂tuˢ = u.momentum_∂tuˢ
                    momentum_p = u.momentum_p
                    momentum_ĝ = u.momentum_ĝ
                    momentum_Fᵥ = u.momentum_Fᵥ
                    T = u.T
                    temperature_turbulence_κ = p.temperature_turbulence_κ
                    nu = p.nu
                    temperature_continuity_C = u.temperature_continuity_C
                    temperature_continuity_F = u.temperature_continuity_F
                    S = u.S
                    salinity_turbulence_κ = p.salinity_turbulence_κ
                    salinity_continuity_C = u.salinity_continuity_C
                    salinity_continuity_F = u.salinity_continuity_F
                    eos_g = p.eos_g
                    eos_α = p.eos_α
                    eos_β = p.eos_β
                    var"0.5" = 0.5
                    var"-1" = -1.0
                    var"-1" = -1.0
                    var"-1.0" = -1.0
                end
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:537 =#
                var"salinity_continuity_•11" = Decapodes.get_tmp(var"__salinity_continuity_•11", u)
                var"momentum_•14" = Decapodes.get_tmp(var"__momentum_•14", u)
                var"momentum_•16" = Decapodes.get_tmp(var"__momentum_•16", u)
                var"salinity_continuity_•9" = Decapodes.get_tmp(var"__salinity_continuity_•9", u)
                var"momentum_•19" = Decapodes.get_tmp(var"__momentum_•19", u)
                var"momentum_•20" = Decapodes.get_tmp(var"__momentum_•20", u)
                var"temperature_continuity_•7" = Decapodes.get_tmp(var"__temperature_continuity_•7", u)
                var"temperature_continuity_•9" = Decapodes.get_tmp(var"__temperature_continuity_•9", u)
                var"temperature_continuity_•11" = Decapodes.get_tmp(var"__temperature_continuity_•11", u)
                var"salinity_continuity_•7" = Decapodes.get_tmp(var"__salinity_continuity_•7", u)
                var"momentum_•18" = Decapodes.get_tmp(var"__momentum_•18", u)
                SD = Decapodes.get_tmp(__SD, u)
                temperature_FD = Decapodes.get_tmp(__temperature_FD, u)
                var"temperature_continuity_•5" = Decapodes.get_tmp(var"__temperature_continuity_•5", u)
                var"temperature_continuity_•4" = Decapodes.get_tmp(var"__temperature_continuity_•4", u)
                var"temperature_continuity_•3" = Decapodes.get_tmp(var"__temperature_continuity_•3", u)
                var"temperature_continuity_•2" = Decapodes.get_tmp(var"__temperature_continuity_•2", u)
                SD = Decapodes.get_tmp(__SD, u)
                salinity_FD = Decapodes.get_tmp(__salinity_FD, u)
                var"salinity_continuity_•5" = Decapodes.get_tmp(var"__salinity_continuity_•5", u)
                var"salinity_continuity_•4" = Decapodes.get_tmp(var"__salinity_continuity_•4", u)
                var"salinity_continuity_•3" = Decapodes.get_tmp(var"__salinity_continuity_•3", u)
                var"salinity_continuity_•2" = Decapodes.get_tmp(var"__salinity_continuity_•2", u)
                var"momentum_•6" = Decapodes.get_tmp(var"__momentum_•6", u)
                temperature_continuity_ċ = Decapodes.get_tmp(__temperature_continuity_ċ, u)
                salinity_continuity_ċ = Decapodes.get_tmp(__salinity_continuity_ċ, u)
                var"momentum_•9" = Decapodes.get_tmp(var"__momentum_•9", u)
                var"momentum_•11" = Decapodes.get_tmp(var"__momentum_•11", u)
                var"momentum_•8" = Decapodes.get_tmp(var"__momentum_•8", u)
                momentum_sum_1 = Decapodes.get_tmp(__momentum_sum_1, u)
                var"momentum_•5" = Decapodes.get_tmp(var"__momentum_•5", u)
                momentum_sum_2 = Decapodes.get_tmp(__momentum_sum_2, u)
                var"momentum_•4" = Decapodes.get_tmp(var"__momentum_•4", u)
                var"momentum_•3" = Decapodes.get_tmp(var"__momentum_•3", u)
                momentum_sum_3 = Decapodes.get_tmp(__momentum_sum_3, u)
                var"momentum_•2" = Decapodes.get_tmp(var"__momentum_•2", u)
                momentum_v̇ = Decapodes.get_tmp(__momentum_v̇, u)
                mul!(var"salinity_continuity_•11", var"GenSim-M_dual_d₀", salinity_continuity_C)
                mul!(var"momentum_•14", var"GenSim-M_dual_d₁", V)
                mul!(var"momentum_•16", var"GenSim-M_dual_d₁", v)
                mul!(var"salinity_continuity_•9", var"GenSim-M_dual_d₀", S)
                mul!(var"momentum_•19", var"GenSim-M_GenSim-ConMat_1", momentum_uˢ)
                mul!(var"momentum_•20", var"GenSim-M_dual_d₀", momentum_p)
                var"temperature_turbulence_•2" = Δᵈ₁(v)
                var"temperature_turbulence_•1" = Δᵈ₀(T)
                mul!(var"temperature_continuity_•7", var"GenSim-M_dual_d₀", T)
                mul!(var"temperature_continuity_•9", var"GenSim-M_dual_d₀", T)
                mul!(var"temperature_continuity_•11", var"GenSim-M_dual_d₀", temperature_continuity_C)
                var"salinity_turbulence_•2" = Δᵈ₁(v)
                var"salinity_turbulence_•1" = Δᵈ₀(S)
                mul!(var"salinity_continuity_•7", var"GenSim-M_dual_d₀", S)
                var"momentum_•7" = ℒ₁(v, v)
                var"momentum_•10" = ι₁₁(v, v)
                var"momentum_•12" = ι₁₁(v, V)
                var"momentum_•13" = ι₁₂(v, var"momentum_•14")
                var"momentum_•15" = ι₁₂(V, var"momentum_•16")
                var"momentum_•18" .= momentum_f .- var"momentum_•19"
                var"momentum_•17" = var"momentum_•18" ∧ᵖᵈ₀₁ v
                SD .= nu .* var"temperature_turbulence_•2"
                temperature_FD .= temperature_turbulence_κ .* var"temperature_turbulence_•1"
                var"temperature_continuity_•6" = ι₁₁(v, var"temperature_continuity_•7")
                var"temperature_continuity_•5" .= var"-1" .* var"temperature_continuity_•6"
                var"temperature_continuity_•8" = ι₁₁(V, var"temperature_continuity_•9")
                var"temperature_continuity_•4" .= var"temperature_continuity_•5" .- var"temperature_continuity_•8"
                var"temperature_continuity_•10" = ι₁₁(v, var"temperature_continuity_•11")
                var"temperature_continuity_•3" .= var"temperature_continuity_•4" .- var"temperature_continuity_•10"
                var"temperature_continuity_•2" .= var"temperature_continuity_•3" .- temperature_FD
                SD .= nu .* var"salinity_turbulence_•2"
                salinity_FD .= salinity_turbulence_κ .* var"salinity_turbulence_•1"
                var"salinity_continuity_•6" = ι₁₁(v, var"salinity_continuity_•7")
                var"salinity_continuity_•5" .= var"-1" .* var"salinity_continuity_•6"
                var"salinity_continuity_•8" = ι₁₁(V, var"salinity_continuity_•9")
                var"salinity_continuity_•4" .= var"salinity_continuity_•5" .- var"salinity_continuity_•8"
                var"salinity_continuity_•10" = ι₁₁(v, var"salinity_continuity_•11")
                var"salinity_continuity_•3" .= var"salinity_continuity_•4" .- var"salinity_continuity_•10"
                var"salinity_continuity_•2" .= var"salinity_continuity_•3" .- salinity_FD
                var"eos_•3" = eos_α .* T
                var"eos_•1" = eos_β .* S
                var"eos_•2" = var"eos_•3" .- var"eos_•1"
                b = eos_g .* var"eos_•2"
                var"momentum_•6" .= var"-1.0" .* var"momentum_•7"
                temperature_continuity_ċ .= (.+)(var"temperature_continuity_•2", temperature_continuity_F)
                salinity_continuity_ċ .= (.+)(var"salinity_continuity_•2", salinity_continuity_F)
                mul!(var"momentum_•9", var"GenSim-M_dual_d₀", var"momentum_•10")
                mul!(var"momentum_•11", var"GenSim-M_dual_d₀", var"momentum_•12")
                var"momentum_•8" .= var"0.5" .* var"momentum_•9"
                var"momentum_•21" = b ∧ᵈᵈ₀₁ momentum_ĝ
                momentum_sum_1 .= (.+)(var"momentum_•6", var"momentum_•8")
                var"momentum_•5" .= momentum_sum_1 .- var"momentum_•11"
                momentum_sum_2 .= (.+)(var"momentum_•5", var"momentum_•13", var"momentum_•15")
                var"momentum_•4" .= momentum_sum_2 .- var"momentum_•17"
                var"momentum_•3" .= var"momentum_•4" .- var"momentum_•20"
                momentum_sum_3 .= (.+)(var"momentum_•3", var"momentum_•21")
                var"momentum_•2" .= momentum_sum_3 .- SD
                momentum_v̇ .= (.+)(var"momentum_•2", momentum_∂tuˢ, momentum_Fᵥ)
                #= /home/runner/work/Decapodes.jl/Decapodes.jl/src/simulation.jl:538 =#
                getproperty(du, :v) .= momentum_v̇
                getproperty(du, :T) .= temperature_continuity_ċ
                getproperty(du, :S) .= salinity_continuity_ċ
            end
    end
end