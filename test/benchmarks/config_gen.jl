using TOML

data = Dict()

fields = ["float_type", "code_target", "resolution"]
stages = ["Setup", "Mesh", "Simulate", "Solve"]

field_data = Dict("fields" => join(fields, ","), "stages" => join(stages, ","))
push!(data, string(0) => field_data)

count = 1
for float in ["Float32", "Float64"]
  for target in ["CPUTarget", "CPUTarget"] # Change one to CUDA
    for resolution in [5, 2, 1]
      datum = Dict(fields[1] => float,
                  fields[2] => target,
                  fields[3] => resolution)
      push!(data, string(count) => datum)
      count += 1
    end
  end
end

open("benchmarks_config_heat.toml", "w") do io
  TOML.print(io, data)
end
