# Benchmarks

## Heat
| Task ID |                   statsfile |                        benchfile | resolution | code_target | float_type | Setup Median time | Mesh Median time | Simulate Median time | Solve Median time |   nf |
|---------|-----------------------------|----------------------------------|------------|-------------|------------|-------------------|------------------|----------------------|-------------------|------|
|       3 |  stats_heat_cpu_test_3.jld2 |  benchmarks_heat_cpu_test_3.json |          1 |   CPUTarget |    Float32 |        0.00440888 |         0.279692 |           0.00309735 |          0.558636 | 9327 |
|       6 |  stats_heat_cpu_test_6.jld2 |  benchmarks_heat_cpu_test_6.json |          1 |   CPUTarget |    Float64 |        0.00447047 |         0.324688 |           0.00329365 |           0.63329 | 9297 |
|       3 | stats_heat_cuda_test_3.jld2 | benchmarks_heat_cuda_test_3.json |          1 |  CUDATarget |    Float32 |         0.0052702 |         0.382094 |           0.00457613 |           0.79329 | 9321 |
|       6 | stats_heat_cuda_test_6.jld2 | benchmarks_heat_cuda_test_6.json |          1 |  CUDATarget |    Float64 |        0.00491387 |         0.386088 |            0.0046531 |          0.801826 | 9297 |
|       2 |  stats_heat_cpu_test_2.jld2 |  benchmarks_heat_cpu_test_2.json |          2 |   CPUTarget |    Float32 |         0.0044489 |        0.0662042 |          0.000771337 |         0.0368417 | 2409 |
|       5 |  stats_heat_cpu_test_5.jld2 |  benchmarks_heat_cpu_test_5.json |          2 |   CPUTarget |    Float64 |        0.00449127 |        0.0769773 |           0.00080348 |         0.0415073 | 2373 |
|       2 | stats_heat_cuda_test_2.jld2 | benchmarks_heat_cuda_test_2.json |          2 |  CUDATarget |    Float32 |        0.00496288 |        0.0869387 |           0.00202888 |          0.205252 | 2409 |
|       5 | stats_heat_cuda_test_5.jld2 | benchmarks_heat_cuda_test_5.json |          2 |  CUDATarget |    Float64 |         0.0049011 |        0.0873349 |           0.00192332 |          0.197789 | 2373 |
|       1 |  stats_heat_cpu_test_1.jld2 |  benchmarks_heat_cpu_test_1.json |          5 |   CPUTarget |    Float32 |        0.00444552 |        0.0101653 |          0.000141233 |        0.00136106 |  471 |
|       4 |  stats_heat_cpu_test_4.jld2 |  benchmarks_heat_cpu_test_4.json |          5 |   CPUTarget |    Float64 |        0.00448736 |        0.0124717 |          0.000148126 |        0.00146659 |  435 |
|       1 | stats_heat_cuda_test_1.jld2 | benchmarks_heat_cuda_test_1.json |          5 |  CUDATarget |    Float32 |        0.00492554 |        0.0136025 |           0.00122809 |         0.0408824 |  471 |
|       4 | stats_heat_cuda_test_4.jld2 | benchmarks_heat_cuda_test_4.json |          5 |  CUDATarget |    Float64 |        0.00524847 |        0.0150585 |           0.00122604 |         0.0378521 |  435 |
