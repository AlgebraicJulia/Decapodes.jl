module PetriNetsTest

  using Catlab
  using Catlab.Present
  using Catlab.Programs
  using CombinatorialSpaces
  using CombinatorialSpaces.ExteriorCalculus
  using GeometryBasics
  using LinearAlgebra
  using AlgebraicPetri
  using Test

  using Decapodes.Simulations
  using Decapodes.Diagrams
  using Decapodes.Schedules
  using Decapodes.Examples
  using Decapodes.OpenDiagrams
  using Decapodes.PetriNets
  using Random
  @present DiffusionSpace2D(FreeExtCalc2D) begin
    X::Space
    k::Hom(Form1(X), Form1(X)) # diffusivity of space, usually constant (scalar multiplication)
    proj₁_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
    proj₂_⁰⁰₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
    sum₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
    prod₀::Hom(Form0(X) ⊗ Form0(X), Form0(X))
  end

  SIRD = LabelledReactionNet{Float64, Float64}([:S=>0.0, :I=>0.0, :R=>0.0, :D=>0.0],
                                               (:inf=>0.01)=>((:S,:I)=>(:I,:I)),
                                               (:rec=>0.001)=>(:I=>:R),
                                               (:death=>0.01)=>(:I=>:D))

  expand_pres!(DiffusionSpace2D, SIRD)

  Diffusion = @decapode DiffusionSpace2D begin
    (C, Ċ)::Form0{X}

    Ċ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(k(d₀{X}(C)))))
  end

  Superposition = @decapode DiffusionSpace2D begin
    (C, Ċ₁, Ċ₂)::Form0{X}

    ∂ₜ{Form0{X}}(C) == Ċ₁ + Ċ₂
  end


  sird_dec = pn2dec(DiffusionSpace2D, SIRD)

  compose_epi = @relation (S,I,R) begin
    epi(S,I,R, D, Ṡ₁,İ₁, Ṙ₁, Ḋ₁)
    diff(S, Ṡ₂)
    diff(I, İ₂)
    diff(R, Ṙ₂)
    diff(D, Ḋ₂)
    superpos(Ṡ₁,Ṡ₂, S)
    superpos(İ₁,İ₂, I)
    superpos(Ṙ₁,Ṙ₂, R)
    superpos(Ḋ₁,Ḋ₂, D)
  end


  composed_epi = oapply(compose_epi, vcat(
    [OpenDiagram(sird_dec, [:S, :I, :R, :D, :Ṡ, :İ, :Ṙ, :Ḋ])],
    [OpenDiagram(Diffusion, [:C, :Ċ]) for i in 1:ns(SIRD)],
    [OpenDiagram(Superposition, [:Ċ₁, :Ċ₂, :C]) for i in 1:ns(SIRD)]
    ));

  res = diag2dwd(composed_epi.functor, in_vars = [:S, :I, :R, :D])

  s = EmbeddedDeltaSet2D{Bool, Point{3,Float64}}()
  points = [(0,0,0),(0,0,1),(0,1,0),(0,1,1)]

  add_vertices!(s, 4, point=points)
  glue_sorted_triangle!(s, 1,2,3)
  glue_sorted_triangle!(s, 2,3,4)

  sd = dual(s);

  funcs = sym2func(sd)
  funcs[:k] = Dict(:operator => 1.0 * I(ne(sd)), :type => MatrixFunc())
  merge!(funcs, gen_functions(SIRD, s))

  Examples.contract_matrices!(res, funcs)

  sim, _ = gen_sim(res, funcs, sd);

  Random.seed!(42)
  u  = rand(nv(s)*4)
  du = zeros(Float64,nv(s)*4)

  tempu = copy(u)
  dt = 0.01
  for i in 1:1000
    sim(du, tempu, [],0)
    tempu .+= du * dt
    # Check that population is conserved
  end
  @test all(tempu .!= u)

  mass_diff = sum(vcat([⋆(Val{0}, sd)*(tempu[(1:nv(s)) .+ (i-1) * nv(s)]) for i in 1:4]...)) -
              sum(vcat([⋆(Val{0}, sd)*(u[(1:nv(s)) .+ (i-1) * nv(s)]) for i in 1:4]...))
  @test mass_diff < 1e-6
end
