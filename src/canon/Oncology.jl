module Oncology 

using DiagrammaticEquations
using DiagrammaticEquations.Deca

# for testing
using Decapodes
using ComponentArrays
using CombinatorialSpaces

using ..Canon 
using Markdown

@docapode("Logistic"
  ,"https://www.google.com"
  ,"Eq. 34 from Yi et al.
  A Review of Mathematical Models for Tumor Dynamics and Treatment Resistance
  Evolution of Solid Tumors,
  with f given as logistic growth. (Eq. 5)"
  ,logistic
  ,begin
    C::Form0
    (Dif, Kd, Cmax)::Constant

    fC == C * (1 - C / Cmax)
    ∂ₜ(C) == Dif * Δ(C) + fC - Kd * C
  end
)

@docapode("Gompertz"
  ,"https://www.google.com"
  ,"Eq. 34 from Yi et al.
  A Review of Mathematical Models for Tumor Dynamics and Treatment Resistance
  Evolution of Solid Tumors,
  with f given as Gompertz growth. (Eq. 6)"
  ,gompertz
  ,begin
  C::Form0
  (Dif, Kd, Cmax)::Constant

  fC == C * ln(Cmax / C)
  ∂ₜ(C) == Dif * Δ(C) + fC - Kd * C
end)

## Another way to account for angiogenesis effect on tumor
## growth is by assuming the carrying capacity of the tumor is
## determined by the effective tumor vascular support that is
## in turn affected by the tumor volume (Eqs. 15 and 16). -Review of...

end
