module DataAggr

using DrWatson
@quickactivate :benchmarks

using Test
using BenchmarkTools

include(helpersdir("data_aggr_helper.jl"))

@testset "Benchmark names" begin
    @test get_benchmark_headername("Solve", "Maximum", "Time") == "Solve Maximum Time"
end

@testset "Adding debug data" begin
    data_row_debug = Dict{String, Any}()
    test_key = "1"
    arch = "cpu"
    debug_data = add_debug_simdata!(data_row_debug, test_key, arch)

    debug_vals = values(debug_data)
    @test test_key in debug_vals
    @test get_statsfile_name(test_key, arch) in debug_vals
    @test get_benchfile_name(test_key, arch) in debug_vals
end

@testset "Adding benchmark data" begin
    test_suite = BenchmarkGroup()
    test_key = "1"
    for solve in get_solver_stages()
        test_suite[test_key][solve] = @benchmarkable rand(100) samples=10
    end
    test_run = run(test_suite)

    data_row_bench = Dict{String, Any}()
    median_test_run = median(test_run[test_key])
    add_trial_data!(data_row_bench, median_test_run, "Median")

    for solve in get_solver_stages()
        for stat in ["Time", "Mem"]
            @test get_benchmark_headername(solve, "Median", stat) in keys(data_row_bench)
        end
    end
end

@testset "Add solver stats data" begin
    pseudo_stats = Dict{String, Any}("steps" => 10, "evals" => 100)
    data_row_stats = Dict{String, Any}()
    add_solver_stats_data!(data_row_stats, pseudo_stats)

    @test data_row_stats["steps"] == 10
    @test data_row_stats["evals"] == 100
end

end