module pathfiles

using DrWatson
@quickactivate "benchmarks"

using Test

include(srcdir("paths.jl"))

@testset "Config name generation" begin
    let name = "heat"; arch = "cpu"
        @test "heat_cpu.toml" == get_configname(name, arch)
    end
end

end