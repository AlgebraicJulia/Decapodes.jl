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

import Decapodes: SummationCall

@testset "Test SummationCall" begin
  # Test equality, 2 inputs
  @test Expr(SummationCall(EQUALS, [:x, :y], :z)) == :(z = (.+)(x, y))

  # Test equality, 3 inputs
  @test Expr(SummationCall(EQUALS, [:x, :y, :w], :z)) == :(z = (.+)(x, y, w))

  # Test equality, 1 input
  @test Expr(SummationCall(EQUALS, [:x], :z)) == :(z = (.+)(x))

  # Test equality, 0 inputs
  @test Expr(SummationCall(EQUALS, [], :z)) == :(z = (.+)())
  
  # Test broadcast equality, 33 inputs
  @test Expr(SummationCall(EQUALS, fill(:x, 33), :z)) == Meta.parse("z = sum([" * foldl(*, fill("x, ", 32)) * "x])")

  # Test broadcast equality, 2 inputs
  @test Expr(SummationCall(DOT_EQUALS, [:x, :y], :z)) == :(z .= (.+)(x, y))

  # Test broadcast equality, 3 inputs
  @test Expr(SummationCall(DOT_EQUALS, [:x, :y, :w], :z)) == :(z .= (.+)(x, y, w))

  # Test broadcast equality, 1 input
  @test Expr(SummationCall(DOT_EQUALS, [:x], :z)) == :(z .= (.+)(x))

  # Test broadcast equality, 0 inputs
  @test Expr(SummationCall(DOT_EQUALS, [], :z)) == :(z .= (.+)())
  
  # Test broadcast equality, 33 inputs
  @test Expr(SummationCall(DOT_EQUALS, fill(:x, 33), :z)) == Meta.parse("z .= sum([" * foldl(*, fill("x, ", 32)) * "x])")
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

    # Test correct data for dimension 1 FixedSizeDiffCache
    for form in [:Form0, :Form1, :DualForm1, :DualForm0]
      test_prealloc_tools(AllocVecCall(:V, form, 1, type, CPUTarget()))
    end

    # Test correct data for dimension 2 FixedSizeDiffCache
    for form in [:Form0, :Form1, :Form2, :DualForm2, :DualForm1, :DualForm0]
      test_prealloc_tools(AllocVecCall(:V, form, 2, type, CPUTarget()))
    end
  end

  for type in [Float32, Float64]

    # Test correct data for dimension 1 CuVector
    for form in [:Form0, :Form1, :DualForm1, :DualForm0]
      test_vector_cache(AllocVecCall(:V, form, 1, type, CUDATarget()))
    end

    # Test correct data for dimension 1 CuVector
    for form in [:Form0, :Form1, :Form2, :DualForm2, :DualForm1, :DualForm0]
      test_vector_cache(AllocVecCall(:V, form, 2, type, CUDATarget()))
    end
  end
end

#####################
# Hooking Code Test #
#####################

# This is just a small test to see that hooking as a general concept works
import Decapodes: hook_AVC_caching # TODO: Remove this import since this should eventually be exported
import Decapodes: AbstractGenerationTarget
struct MYTESTTarget <: CPUBackend end

function hook_AVC_caching(c::AllocVecCall, resolved_form::Symbol, ::MYTESTTarget)
  :(Testing)
end

@testset "Test hook_AVC_caching" begin
  @test Expr(AllocVecCall(:V, :Form0, 1, Float64, MYTESTTarget())) == :Testing
  @test Expr(AllocVecCall(:V, :Form0, 1, Float64, CPUTarget())) != :Testing
  @test Expr(AllocVecCall(:V, :Form0, 1, Float64, CUDATarget())) != :Testing
end

#####################
# compile_env Tests #
#####################

import Decapodes: compile_env, InvalidCodeTargetException

@testset "Test compile_env" begin

  # Test that error throws on unknown code target
  let d = @decapode begin end
    struct BadTarget <: AbstractGenerationTarget end
    @test_throws InvalidCodeTargetException compile_env(d, Set{Symbol}([:test]), Symbol[], BadTarget())
  end

end

#######################
# get_vars_code Tests #
#######################

import Decapodes: get_vars_code, AmbiguousNameException

@testset "Test get_vars_code" begin

  # Test that constants parse correctly
  let d = @decapode begin end
    inputs = [:C]
    add_parts!(d, :Var, 1, name=inputs, type=[:Constant])
    @test get_vars_code(d, inputs, Float64, CPUTarget()).args[begin+1] == :(C = __p__.C)
  end

  # Test that parameters parse correctly
  let d = @decapode begin end
    inputs = [:P]
    add_parts!(d, :Var, 1, name=inputs, type=[:Parameter])
    @test get_vars_code(d, inputs, Float64, CPUTarget()).args[begin+1] == :(P = __p__.P(__t__))
  end

  # TODO: Remove when Literals are not parsed as symbols anymore
  # Test that literals parse correctly
  let d = @decapode begin end
    inputs = [Symbol("2")]
    add_parts!(d, :Var, 1, name=inputs, type=[:Literal])
    @test get_vars_code(d, inputs, Float64, CPUTarget()).args[begin+1] == :(var"2" = 2.0)
  end

  # Test that all forms parse correctly
  for form in [:Form0, :Form1, :Form2, :DualForm0, :DualForm1, :DualForm2]
    let d = @decapode begin end
      inputs = [:F]
      add_parts!(d, :Var, 1, name=inputs, type=[form])
      @test get_vars_code(d, inputs, Float64, CPUTarget()).args[begin+1] == :(F = __u__.F)
    end
  end

  # Test that duplicated names fails
  let d = @decapode begin end
    add_parts!(d, :Var, 2, name=[:A, :A], type=[:Constant, :Constant])
    @test_throws AmbiguousNameException get_vars_code(d, [:A], Float64, CPUTarget())
  end

  # Test invalid input var names fails
  let d = @decapode begin end
    add_parts!(d, :Var, 2, name=[:A, :A], type=[:Constant, :Constant])
    @test_throws AmbiguousNameException get_vars_code(d, [:test], Float64, CPUTarget())
  end

  # Test that duplicated names depends only on names
  let d = @decapode begin end
    add_parts!(d, :Var, 2, name=[:A, :A], type=[:Constant, :Literal])
    @test_throws AmbiguousNameException get_vars_code(d, [:test], Float64, CPUTarget())
  end
end

#######################
# Stub Handling Tests #
#######################

import Decapodes: add_stub, get_stub, InvalidStubException, NO_STUB_RETURN

@testset "Stub handling" begin
  reg_stub = :Test
  empty_stub = Symbol("")
  uni_stub = Symbol("Τϵστ")
  tricky_stub = Symbol("Τest")

  # Test add_stub fuctionality
  let my_stub = :MyStub
    @test :Test_MyStub == add_stub(reg_stub, my_stub)
    @test_throws InvalidStubException add_stub(empty_stub, my_stub)
    @test_throws InvalidStubException add_stub(uni_stub, my_stub)
    @test_throws InvalidStubException add_stub(tricky_stub, my_stub)
    @test_throws InvalidStubException add_stub(NO_STUB_RETURN, my_stub)
  end

  # Test get_stub fuctionality
  let my_stub = :MyStub
    @test reg_stub == get_stub(add_stub(reg_stub, my_stub))
    @test_throws InvalidStubException get_stub(Symbol("_Test"))
    @test_throws InvalidStubException get_stub(Symbol("_"))
    @test NO_STUB_RETURN == get_stub(Symbol(""))
    @test NO_STUB_RETURN == get_stub(Symbol("Test"))
  end

end


########################
# gensim Fuzzing Tests #
########################

import Decapodes: UnsupportedDimensionException, UnsupportedStateeltypeException

@testset "Gensim Fuzzing" begin
  let d = @decapode begin end
    @test_throws UnsupportedDimensionException gensim(d, [:test], dimension = 3, stateeltype = Float64, code_target = CPUTarget())
    @test_throws UnsupportedStateeltypeException gensim(d, [:test], dimension = 2, stateeltype = Int64, code_target = CPUTarget())

    @test_throws UnsupportedDimensionException gensim(d, [:test], dimension = 3, stateeltype = Float64, code_target = CUDATarget())
    @test_throws UnsupportedStateeltypeException gensim(d, [:test], dimension = 2, stateeltype = Int64, code_target = CUDATarget())
  end
end

end
