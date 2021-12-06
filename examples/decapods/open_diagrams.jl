using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphics, Catlab.Programs
using CombinatorialSpaces.ExteriorCalculus
using Catlab.Graphics.Graphviz
using Decapods.OpenDiagrams


draw(d) = to_graphviz(d, box_labels=:name, junction_labels=:variable,
graph_attrs=Dict(:start => "2", :overlap=>"scale"), port_labels=true)

using Catlab.Graphs
using Catlab.Graphs.BasicGraphs
function to_graph(J::FinCat)
    g = BasicGraphs.Graph()
    copy_parts!(g, graph(J))
    return g
end

Graphics.to_graphviz(F::FinFunctor; kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(F; kw...))

function GraphvizGraphs.to_graphviz_property_graph(F::FinFunctor; kw...)
    simplify_vname(s) = begin
      table = Dict(
    "Form0(X)" => "Ω₀",
    "Form1(X)" => "Ω₁",
    "Form2(X)" => "Ω₂",
    "DualForm0(X)" => "Ω̃₀",
    "DualForm1(X)" => "Ω̃₁",
    "DualForm2(X)" => "Ω̃₂",
    "otimes(Form0(X),Form0(X))" => "Ω₀²",
    "otimes(Form1(X),Form1(X))" => "Ω₁²",
    "otimes(Form1(X),DualForm2(X))" => "Ω₁×Ω̃₂",
    "otimes(Form1(X),Form1(X),Form1(X))" => "Ω₁³",
    "otimes(Form1(X),Form1(X),Form1(X),Form1(X))" => "Ω₁⁴" 
    )
    if string(s) in keys(table)
        return table[string(s)]
    else
        println(string(s))
        return string(s)
    end
    end
    
    simplify_ename(s) = begin
        b = IOBuffer()
        show_unicode(b, s)
        return replace(String(take!(b)), r"{X}"=>"")
    end
    
    J = dom(F)
    G = graph(J)
    pg = GraphvizGraphs.to_graphviz_property_graph(to_graph(J); kw...)
    for v in vertices(G)
        lᵥ = G[v, :vname]
        tᵥ = simplify_vname(F.ob_map[v])
        set_vprops!(pg, v, Dict(:label => "$(lᵥ):$tᵥ"))
    end
    for e in edges(G)
        tₑ = F.hom_map[e]
        set_eprops!(pg, e, Dict(:label => "$(simplify_ename(tₑ))"))
    end
    pg
end


Graphics.to_graphviz(cosp::OpenDiagram; kw...) =
to_graphviz(GraphvizGraphs.to_graphviz_property_graph(cosp; kw...))


function GraphvizGraphs.to_graphviz_property_graph(cosp::OpenDiagram; kw...)
    pg = GraphvizGraphs.to_graphviz_property_graph(cosp.functor; kw...)
    label(I, l, i) = length(dom(l)) > 1 ? "$(i):$I" : "$I"
    for (I,l) in enumerate(legs(cosp.cospan))
        for i in dom(l)
            v = add_vertex!(pg)
            set_vprops!(pg, v, Dict(:label=>label(I, l,i), :shape=>"circle", :color=>"blue"))
            e = add_edge!(pg, v, l(i))
            set_eprops!(pg, e, Dict(:style=>"dashed"))
        end
    end
    return pg
end


@present DiffusionSpace2D(FreeExtCalc2D) begin
X::Space
k::Hom(Form1(X), Form1(X)) # diffusivity of space, usually constant (scalar multiplication)
end

FicksLaw2D = @free_diagram DiffusionSpace2D begin
    C::Form0{X}
    dC::Form1{X}
    J::DualForm1{X}
    
    dC == d₀{X}(C)
    J == ⋆₁{X}(k(dC))
end
to_graphviz(OpenDiagram(FicksLaw2D, [:C, :J]), node_labels=true, node_attrs=Dict(:shape=>"oval"), prog="dot")


DiffusionConservation2D = @free_diagram DiffusionSpace2D begin
    (C, Ċ)::Form0{X}
    ϕ::DualForm1{X}
    dϕ::DualForm2{X}
  
    dϕ == dual_d₁{X}(ϕ)
    Ċ == ∂ₜ{Form0{X}}(C)
    Ċ == ⋆₀⁻¹{X}(dϕ)
  end;
  
to_graphviz(OpenDiagram(DiffusionConservation2D, [:C, :Ċ, :ϕ]), node_labels=true, node_attrs=Dict(:shape=>"oval"), prog="dot")
  
compose_diffusion = @relation (C, Ċ, ϕ) begin
  ficks_law(C, ϕ)
  mass_conservation(C, Ċ, ϕ)
end

draw(compose_diffusion)

composed_diffusion = oapply(compose_diffusion, [
    OpenDiagram(FicksLaw2D, [:C, :J]),
    OpenDiagram(DiffusionConservation2D, [:C, :Ċ, :ϕ]),
]);


to_graphviz(composed_diffusion, node_labels=true, node_attrs=Dict(:shape=>"oval"), prog="neato")

@present Flow2DQuantities(FreeExtCalc2D) begin
  X::Space
  k::Hom(Form1(X), Form1(X))    # diffusivity (usually scalar multiplication)
  R₀::Hom(Form0(X), Form0(X))    # Ideal gas constant (usually scalar multiplication)
  kᵥ::Hom(Form1(X), Form1(X))    # viscosity (usually scalar multiplication)
  kₚ::Hom(Form0(X), Form0(X))    # compressibility (usually scalar multiplication)
  m⁻¹::Hom(DualForm1(X), DualForm1(X))    # diffusivity (usually scalar multiplication)
  L₀::Hom(Form1(X)⊗DualForm2(X), DualForm2(X))
  L₁::Hom(Form1(X)⊗DualForm1(X), DualForm1(X))
  L₁′::Hom(Form1(X)⊗Form1(X), Form1(X))
  ∂₀::Hom(Form0(X), Form0(X)) # all_wall boundary condition
  ∂₀₋::Hom(Form0(X), Form0(X)) # left/right wall boundary condition
  ∂₀₊::Hom(Form0(X), Form0(X)) # top/bottom wall boundary condition
  ∂₀ₛ::Hom(Form0(X), Form0(X)) # Sticking boundary condition
  ∂₁ₛ::Hom(Form1(X), Form1(X)) # Sticking boundary condition
  ∂₁ₗ₊::Hom(Form1(X), Form1(X)) # Sticking boundary condition
  ∂₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  dneg₁::Hom(DualForm1(X), DualForm1(X)) # In/Out edge flow boundary condition
  neg₁::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  dneg₀::Hom(DualForm0(X), DualForm0(X)) # In/Out edge flow boundary condition
  neg₀::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  mask₁ₑ::Hom(Form1(X), Form1(X)) # In/Out edge flow boundary condition
  mask₀ₗ::Hom(Form0(X), Form0(X)) # In/Out edge flow boundary condition
  half::Hom(Form0(X), Form0(X)) # half
  i₁′::Hom(Form1(X)⊗Form1(X), Form0(X))
  mult₀::Hom(Form0(X)⊗Form0(X), Form0(X)) # elementwise multiplication
  ∑₀::Hom(Form0(X)⊗Form0(X), Form0(X)) # elementwise addition
  ∑₁⁴::Hom(Form1(X)⊗Form1(X)⊗Form1(X)⊗Form1(X), Form1(X)) # elementwise addition

  proj¹₀::Hom(otimes(Form0(X), Form0(X)), Form0(X))
  proj²₀::Hom(otimes(Form0(X), Form0(X)), Form0(X))
  projᵥ::Hom(otimes(Form1(X), DualForm2(X)), Form1(X))
  projₜ::Hom(otimes(Form1(X), DualForm2(X)), DualForm2(X))
  π¹₁::Hom(Form1(X)⊗Form1(X)⊗Form1(X)⊗Form1(X), Form1(X)) # projections
  π²₁::Hom(Form1(X)⊗Form1(X)⊗Form1(X)⊗Form1(X), Form1(X)) # projections
  π³₁::Hom(Form1(X)⊗Form1(X)⊗Form1(X)⊗Form1(X), Form1(X)) # projections
  π⁴₁::Hom(Form1(X)⊗Form1(X)⊗Form1(X)⊗Form1(X), Form1(X)) # projections
  proj¹₁::Hom(otimes(Form1(X), Form1(X)), Form1(X))
  proj²₁::Hom(otimes(Form1(X), Form1(X)), Form1(X))
  
end

@free_diagram Flow2DQuantities begin
  C::Form0{X}     # concentration
  dC::Form0{X}     # concentration
  T::Form0{X}     # temp
  dT::Form0{X}     # change in temp
  ρ::Form0{X}     # density
  dρ::Form0{X}     # change in density
  V::Form1{X}     # Flow field
  dV::Form1{X}    # Flow field time derivative
  P::DualForm1{X} # Flow momentum
  p::Form0{X}     # pressure field
  dp::Form0{X}     # pressure field time derivative
  ϕ::DualForm1{X} # negative diffusion flux
  bc₀::Form0{X}
  bc₁::Form1{X}
  bc₂::Form1{X}
  bc₃::Form0{X}
end

ficks = @free_diagram Flow2DQuantities begin
  # Ficks law
  T::Form0{X}
  ϕ₁::Form0{X}
#   dual_d₁(X) ⋅ ⋆₀⁻¹(X)
  ϕ₁ == ⋆₀⁻¹{X}(dual_d₁{X}(⋆₁{X}(k(d₀{X}(T)))))
end

to_graphviz(OpenDiagram(ficks, [:T]), node_labels=true)

advection = @free_diagram Flow2DQuantities begin
  # The advection of heat
  T::Form0{X}
  ϕ₂::Form0{X}
  V::Form1{X}
  A::otimes{Form1{X}, DualForm2{X}}
  Tₕ::DualForm2{X}
  projᵥ(A) == V
  projₜ(A) == Tₕ
  Tₕ == ⋆₀{X}(T)
  ϕ₂ == neg₀(⋆₀⁻¹{X}(L₀(A)))
end

to_graphviz(OpenDiagram(advection, [:T, :V, :ϕ₂]), node_labels=true)

super = @free_diagram Flow2DQuantities begin
    # super position of two systems by adding fluxes
    T::Form0{X}
    ϕ₁::Form0{X}
    ϕ₂::Form0{X}
    ϕ₁ϕ₂::otimes{Form0{X}, Form0{X}}
    ϕ::Form0{X}
    ϕ == ∑₀(ϕ₁ϕ₂)
    proj¹₀(ϕ₁ϕ₂) == ϕ₁
    proj²₀(ϕ₁ϕ₂) == ϕ₂
    ϕ == ∂ₜ{Form0{X}}(T)
end

to_graphviz(OpenDiagram(super, [:T, :ϕ₁, :ϕ₂]), node_labels=true)  

thermo = @free_diagram Flow2DQuantities begin  
  ## Equation of State
  T::Form0{X}     # temp
  ρ::Form0{X}     # density
  p::Form0{X}     # pressure field
  ρT::Form0{X}    # product of ρ and T 
  ρₓT::otimes{Form0{X}, Form0{X}}
  proj¹₀(ρₓT) == ρ
  proj²₀(ρₓT) == T
  ρT == mult₀(ρₓT)
  p == R₀(ρT)
end

to_graphviz(OpenDiagram(thermo, [:ρ, :T, :p]), node_labels=true)


navier = @free_diagram Flow2DQuantities begin
    ## Navier Stokes
    ρ::Form0{X}     # density
    ρₕ::DualForm2{X}# hodge of density 
    dρ::Form0{X}     # density derivative
    V::Form1{X}     # Flow field
    V²::otimes{Form1{X}, Form1{X}}     # Flow field
    dV::Form1{X}    # Flow field time derivative total
    dV₁::Form1{X}    # Flow field time derivative component 1
    dV₂::Form1{X}    # Flow field time derivative component 2
    dV₃::Form1{X}    # Flow field time derivative component 3
    dV₄::Form1{X}    # Flow field time derivative component 4
    # P::DualForm1{X} # Flow momentum
    p::Form0{X}     # pressure field
    # dp::Form0{X}     # pressure field time derivative
    B::otimes{Form1{X}, DualForm2{X}}
    projᵥ(B) == V
    ρₕ == ⋆₀{X}(ρ)
    projₜ(B) == ρₕ
    dρ == ∂ₜ{Form0{X}}(ρ) 
    dρ == neg₀(⋆₀⁻¹{X}(L₀(B)))

    dVs::otimes{Form1{X},Form1{X},Form1{X},Form1{X}}
    # dV₁ == (V ⊗ V) ⋅ L₁′ ⋅ neg₁
    dV₁ == neg₁(L₁′(V²))
    dV₂ == kᵥ(Δ₁{X}(V))
    # dV₃ ==  (V ⊗ V) ⋅ i₁′ ⋅ half ⋅ d₀(X)
    dV₃ == d₀{X}(half(i₁′(V²)))
    dV₄ == neg₁(d₀{X}(p))
    dV == ∑₁⁴(dVs)    
    dV == ∂ₜ{Form1{X}}(V)
    π¹₁(dVs) == dV₁
    π²₁(dVs) == dV₂
    π³₁(dVs) == dV₃
    π⁴₁(dVs) == dV₄
    proj¹₁(V²) == V
    proj²₁(V²) == V
end
  
draw_undirected(diagram::OpenDiagram) = to_graphviz(diagram,
  node_labels=true, prog="neato",
  node_attrs=Dict(:shape=>"oval"),
  graph_attrs=Dict(:nodesep=>"0.5", :overlap=>"scale")
)


draw_undirected(OpenDiagram(navier, [:V, :ρ]))

composite_physics = @relation (ρ, p, V, T, ϕ₁, ϕ₂) begin
    ficks(T, ϕ₁)
    advection(T, V, ϕ₂)
    superposition(T, ϕ₁, ϕ₂)
    # diffadv(T, V)
    thermo(ρ, T, p)
    navier(V, ρ, p)
end

draw(composite_physics)
oapply(composite_physics, [
    OpenDiagram(ficks, [:T, :ϕ₁]),
    OpenDiagram(advection, [:T, :V, :ϕ₂]),
    OpenDiagram(super, [:T, :ϕ₁, :ϕ₂]),
    OpenDiagram(thermo, [:ρ, :T, :p]),
    OpenDiagram(navier, [:V, :ρ, :p]),
]) |> draw_undirected

composite_diffadv = @relation (T, V) begin
    ficks(T, ϕ₁)
    advection(T,V, ϕ₂)
    superposition(T, ϕ₁, ϕ₂, ϕ)
end

diffadv = oapply(composite_diffadv, [
  OpenDiagram(ficks, [:T, :ϕ₁]),
  OpenDiagram(advection, [:T, :V, :ϕ₂]),
  OpenDiagram(super, [:T, :ϕ₁, :ϕ₂,:ϕ]),
]) 

diffadv |> draw_undirected



composite_physics = @relation (T,V) begin
    transport(T, V)
    thermo(ρ, T, p)
    fluids(V, ρ, p)
end

draw(composite_physics)
oapply(composite_physics, [
    diffadv,
    OpenDiagram(thermo, [:ρ, :T, :p]),
    OpenDiagram(navier, [:V, :ρ, :p]),
]) |> draw_undirected


draw_undirected(diffadv)
draw_undirected(OpenDiagram(thermo, [:ρ, :T, :p]))
draw_undirected(OpenDiagram(navier, [:V, :ρ, :p]))