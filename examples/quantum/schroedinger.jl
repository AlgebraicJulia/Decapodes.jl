using Catlab, CombinatorialSpaces, Decapodes
using LinearAlgebra, MultiScaleArrays, MLStyle, OrdinaryDiffEq
using GeometryBasics: Point2
Point2D = Point2{Float64}

Schroedinger = @decapode begin
  (i,h,m)::Constant
  V::Parameter
  Psi::Form0
  
  ∂ₜ(Psi) == ((-1 * (h^2)/(2*m))*Δ(Psi) + V * Psi) / (i*h)
end

infer_types!(Schroedinger, op1_inf_rules_1D, op2_inf_rules_1D)
resolve_overloads!(Schroedinger, op1_res_rules_1D, op2_res_rules_1D)

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :Δ => x -> begin
      lap = Δ(0,sd)
      lap * x
    end
    :^ => (x,y) -> x .^ y
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

s′ = EmbeddedDeltaSet1D{Bool, Point2D}()
add_vertices!(s′, 30, point=Point2D.(range(-π/2 + π/32, π/2 - π/32, length=30), 0))
add_edges!(s′, 1:nv(s′)-1, 2:nv(s′))
orient!(s′)
s = EmbeddedDeltaDualComplex1D{Bool, Float64, Point2D}(s′)
subdivide_duals!(s, Circumcenter())

sim = eval(gensim(Schroedinger, dimension=1))
fₘ = sim(s, generate)

Psi = zeros(nv(s))

u₀ = construct(PhysicsState, [VectorForm(Psi)], Float64[], [:Psi])
constants_and_parameters = (
  i = 1, # TODO: Relax assumption of Float64
  V = t -> begin
    z = zeros(nv(s))
    z
  end,
  h = 6.5e-16, # Planck constant in [eV s]
  m = 5.49e-4, # mass of electron in [eV]
)

prob = ODEProblem(fₘ, u₀, (0, 1e-16), constants_and_parameters)
soln = solve(prob, Tsit5())

