module Canon

using DiagrammaticEquations
using DiagrammaticEquations.Deca
using ..Decapodes

function create_pode_expr(t) 
  modelname, source, desc, variable, pode = t
  variable = Symbol(variable)
  Base.remove_linenums!(pode)

  export_stmt = :(export $variable)
  modeldef = sprint.(Base.show_unquoted,pode.args)
  modeldef = join(map(modeldef) do line; "$line
                  " end, "\n")

  docstring = "
  ## $modelname

  [Source]($source)

  $desc

  ### Model

  $modeldef 

  "
  assignment = :($variable = @decapode $pode)
  y = quote 
    $export_stmt
    $assignment
    @doc $docstring $variable
  end
end

macro docapode(name, source, description, variable, pode)
  eval(create_pode_expr((name, source, description, variable, pode)))
end

# @docapode(energy_balance, source, "energy balance equation from Budyko Sellers", energy_balance, 
# begin
#   (Tₛ, ASR, OLR, HT)::Form0
#   (C)::Constant

#   Tₛ̇ == ∂ₜ(Tₛ) 

#   Tₛ̇ == (ASR - OLR + HT) ./ C
# end)


insolation = @decapode begin
  Q::Form0
  cosϕᵖ::Constant

  Q == 450 * cosϕᵖ
end

# warming = @decapode begin
#   (Tₛ)::Form0
#   (A)::Form1

#   A == avg₀₁(5.8282*10^(-0.236 * Tₛ)*1.65e7)
# end


tracer = @decapode begin
  (c,C,F,c_up)::Form0
  (v,V,q)::Form1

  c_up == -1*⋆(L(v,⋆(c))) - ⋆(L(V,⋆(c))) - ⋆(L(v,⋆(C))) - ∘(⋆,d,⋆)(q) + F
end

equation_of_state = @decapode begin
  (b,T,S)::Form0
  (g,α,β)::Constant

  b == g*(α*T - β*S)
end

boundary_conditions = @decapode begin
  (S,T)::Form0
  (Ṡ,T_up)::Form0
  v::Form1
  v_up::Form1
  Ṫ == ∂ₜ(T)
  Ṡ == ∂ₜ(S)
  v̇ == ∂ₜ(v)

  Ṫ == ∂_spatial(T_up)
  v̇ == ∂_noslip(v_up)
end

# Kealy = @decapode begin
# grep  (n,w)::DualForm0
#   dX::Form1
#   (a,ν)::Constant
#   ∂ₜ(w) == a - w - w * n^2 + ν * Δ(w)
# end

Turing_Continuous_Ring = @decapode begin
  (X, Y)::Form0
  (μ, ν, a, b, c, d)::Constant

  ∂ₜ(X) == a*X + b*Y + μ*Δ(X)
  ∂ₜ(Y) == c*X + d*Y + ν*Δ(X)
end

Mohamed_Equation10ForN2 = quote
  𝐮::Form1
  (P, 𝑝ᵈ)::Form0
  (negone, half, μ)::Constant

  ∂ₜ(𝐮) == 𝐮̇

  𝑝ᵈ == P + half * i(𝐮,𝐮)

  𝐮̇ == μ * ∘(d, ⋆, d, ⋆)(𝐮) + (negone)*⋆₁⁻¹(∧₁₀ₚᵈ(𝐮, ⋆(d(𝐮)))) + d(𝑝ᵈ)
end


include("physics.jl")
include("chemistry.jl")
# include("biology.jl")
include("environment.jl")


end # Canon

