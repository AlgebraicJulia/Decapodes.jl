using Krylov, CombinatorialSpaces, LinearAlgebra, Random

function cache_vcycle_lap(s::AbstractVector{<:HasDeltaSet})
  Ls = map(sds) do sd ∇²(0,sd) end
  rs = as_matrix.(fs)./4.0 #4 is the biggest row sum that occurs for triforce, this is not clearly the correct scaling
  ps = transpose.(as_matrix.(fs))
  (Ls, rs, ps)
end

function run_multigrid_vcycles(b,As,rs,ps,steps,cycles=1,alg=cg)
  multigrid_vcycles(zeros(length(b)),b,As,rs,ps,steps,cycles,alg)
end

function vcycle_lap_op(sds::AbstractVector{<:HasDeltaSet})
  vcyc_cache = cache_vcycle_lap(sds)
  x -> run_multigrid_vcycles(x, vcyc_cache..., 3, 10, cg)
end

s = triangulated_grid(1,1,1/4,sqrt(3)/2*1/4,Point3D)
fs = reverse(repeated_subdivisions(1,s,triforce_subdivision_map));
sses = map(fs) do f dom(f).delta_set end
push!(sses,s)
sds = map(sses) do s dualize(s,Circumcenter()) end

our_mesh = first(sds);
lap = ∇²(0,our_mesh);

Random.seed!(1337)
b = lap*rand(nv(our_mesh));

multi_lap = vcycle_lap_op(sds);

u = multi_lap(b);
norm(lap*u-b)/norm(b)

u_test = cg(lap,b,u)[1];
norm(lap*u_test-b)/norm(b)
