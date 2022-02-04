module DecapodesTest

  using Catlab
	using Catlab.Present
  using Catlab.Programs
	using CombinatorialSpaces
	using CombinatorialSpaces.ExteriorCalculus
	using GeometryBasics
  using LinearAlgebra
  using Test

	using Decapodes.Simulations
	using Decapodes.Diagrams
	using Decapodes.Schedules
	using Decapodes.OpenDiagrams
	using Decapodes.Examples
  using Random

	@present DiffusionSpace2D(FreeExtCalc2D) begin
	  X::Space
	  k::Hom(Form1(X), Form1(X)) # diffusivity of space, usually constant (scalar multiplication)
	  proj¹_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
	  proj²_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
	  sum₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
	  prod₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
	end

	Diffusion = @free_diagram DiffusionSpace2D begin
	  (C, Ċ)::Form0{X}

	  Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(k(d₀{X}(C)))))
	end

	Dynamics = @free_diagram DiffusionSpace2D begin
	  (Ċ, C)::Form0{X}

	  Ċ == ∂ₜ{Form0{X}}(C)
	end;

	compose_diff = @relation (C,) begin
	  diff(C, Ċ)
	  dynamics(C, Ċ)
	end

	composed_diff = oapply(compose_diff, [OpenDiagram(Diffusion, [:C, :Ċ]), OpenDiagram(Dynamics, [:C, :Ċ])]);

  res = diag2dwd(composed_diff.functor, in_vars = [:C], out_vars = [:Ċ])

	s = EmbeddedDeltaSet2D{Bool, Point{3,Float64}}()
	points = [(0,0,0),(0,0,1),(0,1,0),(0,1,1)]

	add_vertices!(s, 4, point=points)
	glue_sorted_triangle!(s, 1,2,3)
	glue_sorted_triangle!(s, 2,3,4)

	sd = dual(s);

	funcs = sym2func(sd)
	funcs[:k] = Dict(:operator => 1.0 * I(ne(sd)), :type => MatrixFunc())
	Examples.contract_matrices!(res, funcs)


  sim, _ = gen_sim(res, funcs, sd);

  Random.seed!(42)
  u  = [rand() for i in 1:nv(s)]
  du = zeros(Float64,nv(s))

	tempu = copy(u)
  dt = 0.001
  for i in 10000
    sim(du, tempu, [],0)
    tempu .+= du
    # Check that mass is conserved
  end
  @test all(tempu .!= u)
  @test all(sum(⋆(Val{0}, sd)*tempu)-sum(⋆(Val{0}, sd)*u) < 1e-6)
end
