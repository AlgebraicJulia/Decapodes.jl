module DecapodesTest

  using Catlab
	using Catlab.Present
  using Catlab.Programs
	using CombinatorialSpaces
  using Catlab.WiringDiagrams
	using CombinatorialSpaces.ExteriorCalculus
  using Test

	using Decapodes.Examples

  rules = gen_dec_rules()
  for v in values(rules)
    @test v isa WiringDiagram
  end

end
