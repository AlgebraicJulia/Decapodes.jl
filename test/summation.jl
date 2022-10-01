using Test
using Catlab
using Catlab.Theories
using Catlab.Present
using Catlab.CategoricalAlgebra

using Decapodes
import Decapodes: DecaExpr, eval_eq!


savevizsvg(g, fname::String) = open(fname, "w") do fp
  run_graphviz(fp, to_graphviz(to_graphviz_property_graph(nsdp)), prog="neato", format="svg")
end



NavierStokes = quote
    V::Form1{X}
    V̇::Form1{X}
    G::Form1{X}
    # T::Form0{X}
    ρ::Form0{X}
    ṗ::Form0{X}
    p::Form0{X}
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
parse_decapode(NavierStokes)
nsdp = SummationDecapode(parse_decapode(NavierStokes))
savevizsvg(nsdp, "physics.svg")

NavierStokes_unnested = quote
    V::Form1{X}
    V̇::Form1{X}
    G::Form1{X}
    ρ::Form0{X}
    ṗ::Form0{X}
    p::Form0{X}
    
  
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
parse_decapode(NavierStokes_unnested)
nsdp = SummationDecapode(parse_decapode(NavierStokes_unnested))
savevizsvg(nsdp, "physics_unnested.svg")

