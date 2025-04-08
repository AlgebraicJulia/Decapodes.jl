using ACSets
using CombinatorialSpaces
using BenchmarkTools
using LinearAlgebra
using GeometryBasics: Point2, Point3
using Krylov, KrylovPreconditioners
const Point2D = Point2{Float64}
const Point3D = Point3{Float64}

println("Generating Mesh")
s = triangulated_grid(240, 240, 10, 10, Point2D);
multi_mesh = PrimalGeometricMapSeries(s, binary_subdivision_map, 4);
sd = finest_mesh(multi_mesh);

println("Mesh elements:")
println(" Vertices: $(nv(sd))")
println(" Edges: $(ne(sd))")
println(" Triangles: $(ntriangles(sd))")

println(" Dual vertices: $(nparts(sd, :DualV))")
println(" Dual edges: $(nparts(sd, :DualE))")
println(" Dual triangles: $(nparts(sd, :DualTri))")

lap = ∇²(0, sd)
f = lap * rand(nv(sd));

# Steady state Δu + f = 0 => Δu = -f => u = -∇⁻¹(f)

rel_err(A, x, b; p = 2) = norm(b - A * x, p) / norm(b, p)

# Direct Solver
println("Direct Solver (LU)")
lap_lu = lu(lap)

ds_u = lap_lu \ f;
@btime ds_u = lap_lu \ f;
ds_u = lap_lu \ f;
@show rel_err(lap, ds_u, f)

# Krylov with ilu(0)
println("ILU0 Preconditioned CG")
lap_ilu0 = ilu(lap)

@btime kry_u, stats = cg(lap, f; M = lap_ilu0, ldiv = true);
kry_u, stats = cg(lap, f; M = lap_ilu0, ldiv = true);
@show rel_err(lap, kry_u, f)

cdr2(md::MultigridData) =
  length(md) > 1 ?
    MultigridData(md.operators[1:end-1],md.restrictions[1:end-1],md.prolongations[1:end-1],md.steps[1:end-1]) :
    error("Not enough grids remaining in $md to take the cdr2.")

# CombinatorialSpaces Multigrid
println("Multigrid Solver")
md = MGData(multi_mesh, sd -> ∇²(0, sd), 3, BinarySubdivision());
warm_start = zeros(nv(sd));
@btime mg_u = multigrid_wcycles(warm_start, f, md, 5);
mg_u = multigrid_wcycles(warm_start, f, md, 5);
@show rel_err(lap, mg_u, f)

println("Geometric Multigrid Preconditioned CG")
import LinearAlgebra: ldiv!
function ldiv!(y, A::MultigridData, b)
  y .= multigrid_wcycles(warm_start, b, A, 2);
end

@btime kr_mg_u, kr_mg_stats = cg(lap, f, M = md, ldiv = true);
kr_mg_u, kr_mg_stats = cg(lap, f, M = md, ldiv = true);
@show rel_err(lap, kr_mg_u, f)
