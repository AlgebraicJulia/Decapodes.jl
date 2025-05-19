module Oncology 

using DiagrammaticEquations
using DiagrammaticEquations.Deca

# for testing
using Decapodes
using ComponentArrays
using CombinatorialSpaces

using ..Canon 
using Markdown

@docapode("TumorInvasion"
  ,"https://en.wikipedia.org/wiki/Cancer_cell#Causes"
  ,"[yin_review_2019](@cite), Eq. 35."
  ,invasion
  ,begin
    (C,fC)::Form0
    (Dif, Kd, Cmax)::Constant

    ∂ₜ(C) == Dif * Δ(C) + fC - Kd * C
end)

@docapode("Logistic"
  ,"https://en.wikipedia.org/wiki/Logistic_function"
  ,"[yin_review_2019](@cite), Eq. 5."
  ,logistic
  ,begin
    (C,fC)::Form0
    Cmax::Constant

    fC == C * (1 - C / Cmax)
  end)

@docapode("Gompertz"
  ,"https://en.wikipedia.org/wiki/Gompertz_function"
  ,"[yin_review_2019](@cite), Eq. 6"
  ,gompertz
  ,begin
    (C,fC)::Form0
    Cmax::Constant

    fC == C * ln(Cmax / C)
end)

## Another way to account for angiogenesis effect on tumor
## growth is by assuming the carrying capacity of the tumor is
## determined by the effective tumor vascular support that is
## in turn affected by the tumor volume (Eqs. 15 and 16). -Review of...

end
