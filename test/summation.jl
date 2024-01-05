using Test
using Catlab
using Catlab.Theories
using Catlab.CategoricalAlgebra
using DiagrammaticEquations
using DiagrammaticEquations.Deca

using Decapodes
import Decapodes: DecaExpr, eval_eq!

using Catlab.Graphics
using Catlab.Graphics.Graphviz

NavierStokes = quote
    (V, V̇, G)::Form1{X}
    # T::Form0{X}
    (ρ, ṗ, p)::Form0{X}
    tmp1::Form1{X}
    
    tmp1 == Δ₁(V) + third(d₀(δ₁(V)))
  
    V̇ == neg₁(L₁′(V, V)) + 
          div₁(kᵥ(tmp1), avg₀₁(ρ)) +
          d₀(half(i₁′(V, V))) +
          neg₁(div₁(d₀(p),avg₀₁(ρ))) +
          G
    ∂ₜ(V) == V̇
    ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
    ∂ₜ(p) == ṗ
end

parse_decapode(NavierStokes)
nsdp = SummationDecapode(parse_decapode(NavierStokes))
@test nparts(nsdp, :Σ) == 2
@test nparts(nsdp, :Var) == 26
@test nparts(nsdp, :Op2) == 5
@test nparts(nsdp, :Summand) == 7
to_graphviz_property_graph(nsdp)

NavierStokes_unnested = quote
    (V, V̇, G)::Form1{X}
    (ρ, ṗ, p)::Form0{X}
  
    V̇ == neg₁(L₁′(V, V)) + 
          div₁(kᵥ(Δ₁(V) + third(d₀(δ₁(V)))), avg₀₁(ρ)) +
          d₀(half(i₁′(V, V))) +
          neg₁(div₁(d₀(p),avg₀₁(ρ))) +
          G
    ∂ₜ(V) == V̇
    ṗ == neg₀(⋆₀⁻¹(L₀(V, ⋆₀(p))))# + ⋆₀⁻¹(dual_d₁(⋆₁(kᵨ(d₀(ρ)))))
    ∂ₜ(p) == ṗ
  end

parse_decapode(NavierStokes_unnested)
nsdp = SummationDecapode(parse_decapode(NavierStokes_unnested))
@test nparts(nsdp, :Σ) == 2
@test nparts(nsdp, :Var) == 26
@test nparts(nsdp, :Op2) == 5
@test nparts(nsdp, :Summand) == 7

to_graphviz_property_graph(nsdp)
