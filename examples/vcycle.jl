using ACSets
using CombinatorialSpaces
using ComponentArrays
using Krylov
using LinearAlgebra
using Random
using DiagrammaticEquations
using Decapodes
using MLStyle
using OrdinaryDiffEq

import Decapodes: default_dec_matrix_generate

s = triangulated_grid(1,1,1/4,sqrt(3)/2*1/4,Point3D)
fs = reverse(repeated_subdivisions(3,s,triforce_subdivision_map));

our_mesh = dualize(dom(first(fs)).delta_set, Circumcenter());
lap = ∇²(0,our_mesh);

Random.seed!(1337)
b = lap*rand(nv(our_mesh));

multi_lap = vcycle_lap_op(fs);

u = multi_lap(b);
norm(lap*u-b)/norm(b)

u_test = cg(lap,b,u)[1];
norm(lap*u_test-b)/norm(b)

inv_lap = @decapode begin
  U::Form0
  C::Form1
  C == d(Δ₀⁻¹(U))
end

function generate(fs, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    _ => default_dec_matrix_generate(fs, my_symbol, hodge)
  end
end

sim = eval(gensim(inv_lap; multigrid=true))

f = sim(fs, generate);
u₀ = ComponentArray(U=b)

prob = ODEProblem(f,u₀,(0,1),());
soln = solve(prob, Tsit5(), adaptive=false, dt=1);

soln.u
