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
  ,""
  ,logistic
  ,begin
    r::Constant # intrinsic growth rate
    K::Constant # carrying capacity
    V::Form0

    ∂ₜ(V) == r * V * (1 - V / K)
  end
)

logistic = @decapode begin
  K::Constant
  r::Constant
  V::Form0

  ∂ₜ(V) == r * V * (1 - V / K)
end

logadv = @decapode begin
  K₀::Constant # carrying-capacity
  c::Constant  # waste metabolite process
  r::Constant  # TODO should be contingent on nutrient influx
  K::Form0     # current carrying-capacity
  W::Form0     # waste
  V::Form0
  ∂ₜ(W) == c * V
  K == K₀- W
  ∂ₜ(V) == r * V * (1 - V / K)
end

## Another way to account for angiogenesis effect on tumor
## growth is by assuming the carrying capacity of the tumor is
## determined by the effective tumor vascular support that is
## in turn affected by the tumor volume (Eqs. 15 and 16). -Review of...

# gompertz = @decapode begin
#   {a,b,c}::Constant
#   V::Form0
#   ∂ₜ(V) == a * V * log(b / (V + c))
# end

space = Icosphere(1) |> loadmesh

sim = evalsim(logistic)

f = sim(space, default_dec_generate)

forms = ComponentArray(V = rand(nv(space)))

copied = copy(forms)

params = ComponentArray(r = 1, K = 2)

f(copied, forms, params, 0)


end
