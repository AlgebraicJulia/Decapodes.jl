module TontiTest
  using Decapods.TontiDiagrams
  using CombinatorialSpaces
  using GeometryBasics
  using Random
  using Test
  td = TontiDiagram();
  add_variables!(td, (:C,0,true),  (:dC,0,true),  (:ΔC,1,true),
                 (:∑ϕ,2,false), (:ϕ,1,false))

  s = EmbeddedDeltaSet2D{Bool, Point{3,Float64}}()
  points = [(0,0,0),(0,0,1),(0,1,0),(0,1,1)]

  add_vertices!(s, 4, point=points)
  glue_sorted_triangle!(s, 1,2,3)
  glue_sorted_triangle!(s, 2,3,4)

  sp = Space(s)

  @test gen_form(s, x->1) == [1,1,1,1]

  add_derivatives!(td, sp, :C=>:ΔC, :ϕ =>:∑ϕ)
  add_time_dep!(td, :dC, :C)

  const k = 0.003

  add_transition!(td, [:ΔC], (ϕ,ΔC)->(ϕ .= sp.hodge[1,2] * (k .* ΔC)), [:ϕ])
  add_transition!(td, [:∑ϕ], (dC, ∑ϕ)->(dC .= sp.hodge[2,3] * ∑ϕ), [:dC]);

  data_range, sim= vectorfield(td, sp)

  @test data_range[:C] == (1,4)

  Random.seed!(42)
  u  = gen_form(s,x->rand())
  du = zeros(Float64,nv(s))

  tempu = copy(u)
  dt = 0.001
  for i in 10000
    sim(du, tempu, [],0)
    tempu .+= du
    # Check that mass is conserved
  end
  @test all(tempu .!= u)
  @test all(sum(sp.hodge[1,1]*tempu)-sum(sp.hodge[1,1]*u) < 1e-6)
end
