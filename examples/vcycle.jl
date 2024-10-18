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
  C::Form0
  C == Δ₀⁻¹(U)
end

function generate(fs, my_symbol; hodge=DiagonalHodge())
  op = @match my_symbol begin
    :Δ₀⁻¹ => vcycle_lap_op(fs)
    _ => default_dec_matrix_generate(fs, my_symbol, hodge)
  end
end

# sim = eval(gensim(inv_lap))

answer = []

sim = (meshs, operators, hodge = GeometricHodge())->begin
            begin
              Δ₀⁻¹ = default_dec_matrix_generate(meshs, :Δ₀⁻¹, hodge)
            end
            f(du, u, p, t) = begin
                    begin
                        U = u.U
                    end
                    C = Δ₀⁻¹(U)
                    push!(answer, C)
                    return nothing
                end
        end
f = sim(fs, generate);
u₀ = ComponentArray(U=b)

prob = ODEProblem(f,u₀,(0,1),());
soln = solve(prob, Tsit5(), adaptive=false, dt=1);

soln.u

u_deca = last(answer)
norm(lap*u_deca-b)/norm(b)
norm(lap*u_test-b)/norm(b)
