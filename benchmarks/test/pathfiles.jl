module PathFiles

using DrWatson
@quickactivate :benchmarks

using Test

@testset "File names" begin
  @test "heat_cpu.toml" == get_configname("heat", "cpu")
  @test occursin("stats_1_cpu.jld2", get_statsfile("1", "heat", "cpu"))
  @test occursin("benchmarks_1_cpu.json", get_benchfile("1", "heat", "cpu"))
end

end