using TOML

data = Dict()

count = 1
for float in ["Float32", "Float64"]
  for target in ["CPUTarget", "CPUTarget"] # Change one to CUDA
    datum = Dict("float_type" => float, "code_target" => target)
    push!(data, string(count) => datum)
    count += 1
  end
end

open("benchmarks_config_heat.toml", "w") do io
  TOML.print(io, data)
end
