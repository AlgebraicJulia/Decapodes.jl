using TOML
using MLStyle

function generate_config(sim_name, target)

  code_target = @match target begin
    "cpu" => "CPUTarget"
    "cuda" => "CUDATarget"
    _ => error("Please supply target as either 'cpu' or 'cuda'")
  end

  data = Dict()

  fields = ["float_type", "code_target", "resolution"]
  stages = ["Setup", "Mesh", "Simulate", "Solve"]

  idx = 0
  for float in ["Float32", "Float64"]
    for resolution in [5, 2, 1]
      idx += 1
      datum = Dict(fields[1] => float,
                  fields[2] => code_target,
                  fields[3] => resolution)
      push!(data, string(idx) => datum)
    end
  end

  meta_data = Dict("fields" => join(fields, ","),
                   "stages" => join(stages, ","),
                   "count" => idx)

  open("count_$(target).txt", "w") do file
    write(file, string(idx))
  end

  push!(data, string(0) => meta_data)

  open("$(sim_name)_$(target).toml", "w") do io
    TOML.print(io, data)
  end
end

function update_tracker(file_name)
  idx = isfile(file_name) ? parse(Int64, readline(file_name)) : 0
  open(file_name, "w") do file
    write(file, string(idx + 1))
  end
end

cd(ARGS[1])
generate_config(ARGS[1], ARGS[2])
update_tracker("config_tracker.txt")
