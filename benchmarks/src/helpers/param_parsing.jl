# Collection of functions meant to convert TOML data into a sim-usable format
using MLStyle
using Decapodes

export parse_float_type, parse_code_target, arch_to_code_target

function parse_float_type(float_data)
    @match float_data begin
        "Float32" => Float32
        "Float64" => Float64
        "ComplexF32" => ComplexF32
        "ComplexF64" => ComplexF64
        _ => error("Float data $(float_data) is not valid")
    end
end

function parse_code_target(code_target_data)
    @match code_target_data begin
        "CPUTarget" => CPUTarget()
        "CUDATarget" => CUDATarget()
        _ => error("Code target data $(code_target_data) is not in list [\"CPUTarget\", \"CUDATarget\"]")
    end
end
  
function arch_to_code_target(architecture)
    @match architecture begin
        "cpu" => "CPUTarget"
        "cuda" => "CUDATarget"
        _ => error("Architecture $(architecture) is not in list [\"cpu\", \"cuda\"]")
    end
end
