using Krylov, CombinatorialSpaces, LinearAlgebra, Random

function cache_vcycle_lap(s::HasDeltaSet)
  fs = reverse(repeated_subdivisions(1,s,triforce_subdivision_map));
  sses = map(fs) do f dom(f).delta_set end
  push!(sses,s)
  sds = map(sses) do s dualize(s,Circumcenter()) end
  Ls = map(sds) do sd ∇²(0,sd) end
  rs = as_matrix.(fs)./4.0 #4 is the biggest row sum that occurs for triforce, this is not clearly the correct scaling
  ps = transpose.(as_matrix.(fs))
  (Ls, rs, ps)
end

function run_multigrid_vcycles(b,As,rs,ps,steps,cycles=1,alg=cg)
  multigrid_vcycles(zeros(length(b)),b,As,rs,ps,steps,cycles,alg)
end

function vcycle_lap_op(s::HasDeltaSet)
  vcyc_cache = cache_vcycle_lap(s)
  x -> run_multigrid_vcycles(x, vcyc_cache..., 3, 10, cg)
end

s = triangulated_grid(1,1,1/4,sqrt(3)/2*1/4,Point3D)
sd = dualize(s,Circumcenter())
lap = ∇²(0,sd)

Random.seed!(1337)
b = Ls[1]*rand(nv(sds[1]))

multi_lap = vcycle_lap_op(s)

u = multi_lap(b)
norm(Ls[1]*u-b)/norm(b)
