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
    sim_name = "heat"
    arch = "cpu"
    tag = "testing"
    test_key = "1"
    test_namedata = SimNameData(sim_name, arch, tag, test_key)

    debug_data = add_debug_simdata!(data_row_debug, test_namedata)

    debug_vals = values(debug_data)
    @test test_key in debug_vals
    @test arch in debug_vals
    @test sim_name in debug_vals
    @test statsfile_name(test_namedata) in debug_vals
    @test benchfile_name(test_namedata) in debug_vals
end

@testset "Adding benchmark data" begin
    test_suite = BenchmarkGroup()
    test_key = "1"
    for stage in solver_stages()
        test_suite[test_key][stage] = @benchmarkable rand(100) samples=10
    end
    test_run = run(test_suite)

    data_row_bench = Dict{String, Any}()
    median_test_run = median(test_run[test_key])
    add_trial_data!(data_row_bench, median_test_run, "Median")

    for stage in solver_stages()
        for stat in ["Time", "Mem"]
            @test get_benchmark_headername(stage, "Median", stat) in keys(data_row_bench)
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