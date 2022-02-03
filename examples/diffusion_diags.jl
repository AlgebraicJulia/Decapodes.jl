using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphics, Catlab.Programs
using CombinatorialSpaces.ExteriorCalculus

@present DiffusionSpace2D(FreeExtCalc2D) begin
  X::Space
  k::Hom(Form1(X), Form1(X)) # diffusivity of space, usually constant (scalar multiplication)
  L::Hom(Form0(X)⊗Form0(X), Form0(X))
end

FicksLaw2D = @free_diagram DiffusionSpace2D begin
  C::Form0{X}
  dC::Form1{X}
  J::DualForm1{X}

  dC == d₀{X}(C)
  J == ⋆₁{X}(k(dC))
end

DiffusionConservation2D = @free_diagram DiffusionSpace2D begin
  (C, Ċ)::Form0{X}
  ϕ::DualForm1{X}
  dϕ::DualForm2{X}

  dϕ == dual_d₁{X}(ϕ)
  Ċ == ∂ₜ{Form0{X}}(C)
  Ċ == ⋆₀⁻¹{X}(dϕ)
end;

using Decapodes

compose_diffusion = @relation (C, Ċ, ϕ) begin
  ficks_law(C, ϕ)
  mass_conservation(C, Ċ, ϕ)
end

Decapodes.OpenDiagrams.draw(compose_diffusion)

using Decapodes.OpenDiagrams

composed_diffusion = oapply(compose_diffusion, [
    OpenDiagram(FicksLaw2D, [:C, :J]),
    OpenDiagram(DiffusionConservation2D, [:C, :Ċ, :ϕ]),
]);

res = diag2dwd(composed_diffusion.functor)