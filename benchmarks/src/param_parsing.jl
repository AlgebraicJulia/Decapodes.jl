
function parse_float_type(float_data)
    @match float_data begin
        "Float32" => Float32
        "Float64" => Float64
        "ComplexF32" => ComplexF32
        "ComplexF64" => ComplexF64
        _ => error("Float data $(float_data) is not valid, exiting early") 
    end
end

function parse_code_target(code_target_data)
    @match code_target_data begin
        "CPUTarget" => CPUTarget()
        "CUDATarget" => CUDATarget()
        _ => error("Warning: Code target data $(code_target_data) is not valid, exiting early")
    end
end
  