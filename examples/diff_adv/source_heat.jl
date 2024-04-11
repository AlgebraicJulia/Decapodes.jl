using Catlab
using CombinatorialSpaces
using Decapodes
using DiagrammaticEquations
using Distributions
using MLStyle
using OrdinaryDiffEq
using LinearAlgebra
using ComponentArrays
using CairoMakie
using CUDA
using CUDA.CUSPARSE
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

function show_heatmap(Cdata)
  heatmap(reshape(Cdata, (floor(Int64, sqrt(length(Cdata))), floor(Int64, sqrt(length(Cdata))))))
end

s = triangulated_grid(100,100,1,1,Point3D);
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s);
subdivide_duals!(sd, Circumcenter());

function generate(sd, symbol, hodge=DiagonalHodge())
  op = @match symbol begin
    _ => error("Unmatched operator $my_symbol")
  end
  return op
end

Heat = @decapode begin
  C::Form0
  dX::Form1
  (Dif, c)::Constant
  S::Parameter
  ∂ₜ(C) == Dif * Δ(C) + c * S + ⋆(d(⋆(C ∧ dX)))
end

sim = eval(gensim(Heat, code_target = cuda()))

function simulate(mesh, operators, hodge = GeometricHodge())
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:577 =#
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:578 =#
  begin
      #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:174 =#
      (var"GenSim-M_⋆₀⁻¹", ⋆₀⁻¹) = default_dec_cu_matrix_generate(mesh, :⋆₀⁻¹, hodge)
      (var"GenSim-M_⋆₁", ⋆₁) = default_dec_cu_matrix_generate(mesh, :⋆₁, hodge)
      (var"GenSim-M_dual_d₁", dual_d₁) = default_dec_cu_matrix_generate(mesh, :dual_d₁, hodge)
      (var"GenSim-M_d₀", d₀) = default_dec_cu_matrix_generate(mesh, :d₀, hodge)
      (var"GenSim-M_∧₀₁", ∧₀₁) = default_dec_cu_matrix_generate(mesh, :∧₀₁, hodge)
  end
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:579 =#
  begin
      #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:484 =#
      var"GenSim-M_GenSim-ConMat_1" = var"GenSim-M_⋆₀⁻¹" * var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁" * var"GenSim-M_d₀"
      var"GenSim-ConMat_1" = (x->var"GenSim-M_GenSim-ConMat_1" * x)
      var"GenSim-M_GenSim-ConMat_2" = var"GenSim-M_⋆₀⁻¹" * (var"GenSim-M_dual_d₁" * var"GenSim-M_⋆₁")
      var"GenSim-ConMat_2" = (x->var"GenSim-M_GenSim-ConMat_2" * x)
  end
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:580 =#
  begin
      #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:226 =#
      var"•3" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•2" = CuVector{Float64}(undef, nparts(mesh, :V))
      var"•8" = CuVector{Float64}(undef, nparts(mesh, :E))
      var"•5" = CuVector{Float64}(undef, nparts(mesh, :V))
      Ċ = CuVector{Float64}(undef, nparts(mesh, :V))
  end
  #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:581 =#
  f(du, u, p, t) = begin
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:581 =#
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:582 =#
          begin
              #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:251 =#
              C = u.C
              dX = u.dX
              Dif = p.Dif
              c = p.c
              S = p.S(t)
          end
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:583 =#
          mul!(var"•3", var"GenSim-M_GenSim-ConMat_1", C)
          var"•2" .= Dif .* var"•3"
          var"•4" = c .* S
          var"GenSim-M_∧₀₁"(var"•8", C, dX)
          mul!(var"•5", var"GenSim-M_GenSim-ConMat_2", var"•8")
          Ċ .= (.+)(var"•2", var"•4", var"•5")
          #= c:\Users\georger\Documents\GitHub\Decapodes.jl\src\simulation.jl:584 =#
          getproperty(du, :C) .= Ċ
      end
end
fₘ = simulate(sd, generate, DiagonalHodge())

t_dist = MvNormal([75, 75], 1)
S = CuVector([pdf(t_dist, [p[1], p[2]]) for p in sd[:point]])

constants_and_parameters = (
  c = 5000 * 1000,
  Dif = 20,
  S = t -> S * t
)

C = CUDA.zeros(Float64, nv(sd))
dX = CuVector(sd[:length])

u₀ = ComponentArray(C=C, dX=dX)
tₑ = 15.0

prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5(), save_idxs=[:C])

begin
  frames = 100
  fig = Figure()
  ax = CairoMakie.Axis(fig[1,1])
  msh = CairoMakie.mesh!(ax, s, color=Array(soln(0).C), colormap=:jet, colorrange=extrema(Array(soln(tₑ).C)))
  Colorbar(fig[1,2], msh)
  CairoMakie.record(fig, "Ocillating Heat.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
      msh.color = Array(soln(t).C)
  end
end
