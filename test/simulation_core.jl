module SimulationCoreTest

using ACSets
using Decapodes
using DiagrammaticEquations
using MLStyle
using Test

const EQUALS = :(=)
const DOT_EQUALS =:.=

######################
# Calling Code Tests #
######################

import Decapodes: UnaryCall, add_inplace_stub

@testset "Test UnaryCall" begin
  # Test equality, basic operator
  @test Expr(UnaryCall(:F, EQUALS, :x, :y)) ==  :(y = F(x))

  # Test equality, inplace negation
  @test Expr(UnaryCall(:.-, EQUALS, :x, :y)) ==  :(y = .-x)

  # Test equality, inplacable op
  let inplacable_op = :⋆₁⁻¹
    @test Expr(UnaryCall(inplacable_op, EQUALS, :x, :y)) ==  :(y = (⋆₁⁻¹)(x))
  end

  # Test broadcast equality, basic operator
  @test Expr(UnaryCall(:F, DOT_EQUALS, :x, :y)) == :(mul!(y, F, x))

  # Test broadcast equality, inplace negation
  @test Expr(UnaryCall(:.-, DOT_EQUALS, :x, :y)) == :(y .= .-x)

  # Test broadcast equality, inplace non-matrix method
  let inplace_op = add_inplace_stub(:⋆₁⁻¹)
    @test Expr(UnaryCall(inplace_op, DOT_EQUALS, :x, :y)) == :($inplace_op(y, x))
  end
end

import Decapodes: BinaryCall

@testset "Test BinaryCall" begin
  # Test equality, basic operator
  @test Expr(BinaryCall(:F, EQUALS, :x, :y, :z)) == :(z = F(x, y))

  # Test broadcast equality, inplace operator
  let inplace_operator = add_inplace_stub(:F)
    @test Expr(BinaryCall(inplace_operator, DOT_EQUALS, :x, :y, :z)) == :($inplace_operator(z, x, y))
  end
end

import Decapodes: VarargsCall

@testset "Test VarargsCall" begin
  # Test equality, 2 inputs
  @test Expr(VarargsCall(:F, EQUALS, [:x, :y], :z)) == :(z = F(x, y))

  # Test equality, 3 inputs
  @test Expr(VarargsCall(:F, EQUALS, [:x, :y, :w], :z)) == :(z = F(x, y, w))

  # Test equality, 1 input
  @test Expr(VarargsCall(:F, EQUALS, [:x], :z)) == :(z = F(x))

  # Test equality, 0 inputs
  @test Expr(VarargsCall(:F, EQUALS, [], :z)) == :(z = F())

  # Test broadcast equality, 2 inputs
  @test Expr(VarargsCall(:F, DOT_EQUALS, [:x, :y], :z)) == :(z .= F(x, y))

  # Test broadcast equality, 3 inputs
  @test Expr(VarargsCall(:F, DOT_EQUALS, [:x, :y, :w], :z)) == :(z .= F(x, y, w))

  # Test broadcast equality, 1 input
  @test Expr(VarargsCall(:F, DOT_EQUALS, [:x], :z)) == :(z .= F(x))

  # Test broadcast equality, 0 inputs
  @test Expr(VarargsCall(:F, DOT_EQUALS, [], :z)) == :(z .= F())
end

#######################
# AllocVec Code Tests #
#######################

import Decapodes: AllocVecCall

@testset "Test AllocVecCall" begin

  function convert_form_to_simplex(form::Symbol, dimension::Int)
    @match (form, dimension) begin
      (:Form0, 2) => :V
      (:Form1, 2) => :E
      (:Form2, 2) => :Tri
      (:DualForm0, 2) => :Tri
      (:DualForm1, 2) => :E
      (:DualForm2, 2) => :V

      (:Form0, 1) => :V
      (:Form1, 1) => :E
      (:DualForm0, 1) => :E
      (:DualForm1, 1) => :V
      _ => @error("Unknown form type for form $(form) in dim $(dimension)")
    end
  end

  """
    test_prealloc_tools(alloc_vec::AllocVecCall)

  This function tests that the PreallocationTools usage is correct. Exactly, this tests that
  the form to simplex conversion is correct, that the correct collection is used and that the
  name of the cache is different from the variable. This only works for `FixedSizeDiffCache`.
  """
  function test_prealloc_tools(alloc_vec::AllocVecCall)
    expr = Expr(alloc_vec)
    name = expr.args[begin]
    body = expr.args[end]

    cache = body.args[begin]
    vec = body.args[end]

    simplex_type = vec.args[end].args[end].value
    collection_type = vec.args[begin]

    type_result = convert_form_to_simplex(alloc_vec.form, alloc_vec.dimension)

    @test alloc_vec.name != name

    @test cache == :(Decapodes.FixedSizeDiffCache)

    @test type_result == simplex_type

    if alloc_vec.code_target isa CPUBackend
      @test collection_type == :(Vector{$(alloc_vec.T)})
    elseif alloc_vec.code_target isa CUDABackend
      @test collection_type == :(CuVector{$(alloc_vec.T)})
    else
      @error("Unknown code target: $(alloc_vec.code_target)")
    end
  end

  """
    test_vector_cache(alloc_vec::AllocVecCall)

  This function tests that the vector caching usage is correct. Exactly, this tests that
  the form to simplex conversion is correct, that the correct collection is used and that the
  name of the cache is the same the variable.

  """
  function test_vector_cache(alloc_vec::AllocVecCall)
    expr = Expr(alloc_vec)
    name = expr.args[begin]
    vec = expr.args[end]

    simplex_type = vec.args[end].args[end].value
    collection_type = vec.args[begin]

    type_result = convert_form_to_simplex(alloc_vec.form, alloc_vec.dimension)

    @test alloc_vec.name == name

    @test type_result == simplex_type

    if alloc_vec.code_target isa CPUBackend
      @test collection_type == :(Vector{$(alloc_vec.T)})
    elseif alloc_vec.code_target isa CUDABackend
      @test collection_type == :(CuVector{$(alloc_vec.T)})
    else
      @error("Unknown code target: $(alloc_vec.code_target)")
    end
  end

  for type in [Float32, Float64]
    for form in [:Form0, :Form1, :DualForm1, :DualForm0]
      test_prealloc_tools(AllocVecCall(:V, form, 1, type, CPUTarget()))
    end

    for form in [:Form0, :Form1, :Form2, :DualForm2, :DualForm1, :DualForm0]
      test_prealloc_tools(AllocVecCall(:V, form, 2, type, CPUTarget()))
    end
  end

  for type in [Float32, Float64]
    for form in [:Form0, :Form1, :DualForm1, :DualForm0]
      test_vector_cache(AllocVecCall(:V, form, 1, type, CUDATarget()))
    end

    for form in [:Form0, :Form1, :Form2, :DualForm2, :DualForm1, :DualForm0]
      test_vector_cache(AllocVecCall(:V, form, 2, type, CUDATarget()))
    end
  end
end

#####################
# Hooking Code Test #
#####################

import Decapodes: hook_AVC_caching # TODO: Remove this import since this should be exported

struct MYTESTTarget <: CPUBackend end

function hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::MYTESTTarget)
  :(Testing)
end

@testset "Test hook_AVC_caching" begin
  @test Expr(AllocVecCall(:V, :Form0, 1, Float64, MYTESTTarget())) == :Testing
  @test Expr(AllocVecCall(:V, :Form0, 1, Float64, CPUTarget())) != :Testing
  @test Expr(AllocVecCall(:V, :Form0, 1, Float64, CUDATarget())) != :Testing

end

end
