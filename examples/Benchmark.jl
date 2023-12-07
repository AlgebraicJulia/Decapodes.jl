using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using LinearAlgebra
using BenchmarkTools
using Random
import Decapodes: dec_differential, dec_dual_derivative, dec_boundary, dec_wedge_product, dec_hodge_star, dec_inv_hodge
using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

begin 
    primal_earth = EmbeddedDeltaSet2D("Ico7.obj");
    orient!(primal_earth);
    earth = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(primal_earth);
    subdivide_duals!(earth, Barycenter());
end

begin 
    println("Mesh: " * "Blender Ico Sphere, 7 Subdivisions")
    println("")

    println("Number of primal vertices: ", nv(earth))
    println("Number of primal edges: ", ne(earth))
    println("Number of primal triangles: ", ntriangles(earth))
    println("")

    println("Number of dual vertices: ", nparts(earth, :DualV))
    println("Number of dual edges: ", nparts(earth, :DualE))
    println("Number of dual triangles: ", nparts(earth, :DualTri))
    println("----------------------------------------------------------------")
end

suite = BenchmarkGroup()

begin
    suite["Exterior Derivative"] = BenchmarkGroup()

    suite["Exterior Derivative"]["New Form-0"] = @benchmarkable dec_differential(0, $earth)
    suite["Exterior Derivative"]["Old Form-0"] = @benchmarkable d(0, $earth)

    suite["Exterior Derivative"]["New Form-1"] = @benchmarkable dec_differential(1, $earth)
    suite["Exterior Derivative"]["Old Form-1"] = @benchmarkable d(1, $earth)
end

begin
    suite["Boundary"] = BenchmarkGroup()

    suite["Boundary"]["New Form-1"] = @benchmarkable dec_boundary(1, $earth)
    suite["Boundary"]["Old Form-1"] = @benchmarkable ∂(1, $earth)

    suite["Boundary"]["New Form-2"] = @benchmarkable dec_boundary(2, $earth)
    suite["Boundary"]["Old Form-2"] = @benchmarkable ∂(2, $earth)
end

begin
    suite["Dual Derivative"] = BenchmarkGroup()

    suite["Dual Derivative"]["New Dual-Form-0"] = @benchmarkable dec_dual_derivative(0, $earth)
    suite["Dual Derivative"]["Old Dual-Form-0"] = @benchmarkable dual_derivative(0, $earth)

    suite["Dual Derivative"]["New Dual-Form-1"] = @benchmarkable dec_dual_derivative(1, $earth)
    suite["Dual Derivative"]["Old Dual-Form-1"] = @benchmarkable dual_derivative(1, $earth)
end

begin
    suite["Diagonal Hodge"] = BenchmarkGroup()

    suite["Diagonal Hodge"]["New Form-0"] = @benchmarkable dec_hodge_star(0, $earth, DiagonalHodge())
    suite["Diagonal Hodge"]["Old Form-0"] = @benchmarkable hodge_star(0, $earth, DiagonalHodge())

    suite["Diagonal Hodge"]["New Form-1"] = @benchmarkable dec_hodge_star(1, $earth, DiagonalHodge())
    suite["Diagonal Hodge"]["Old Form-1"] = @benchmarkable hodge_star(1, $earth, DiagonalHodge())

    suite["Diagonal Hodge"]["New Form-2"] = @benchmarkable dec_hodge_star(2, $earth, DiagonalHodge())
    suite["Diagonal Hodge"]["Old Form-2"] = @benchmarkable hodge_star(2, $earth, DiagonalHodge())

end

begin
    suite["Geometric Hodge"] = BenchmarkGroup()

    suite["Geometric Hodge"]["New Form-1"] = @benchmarkable dec_hodge_star(1, $earth, GeometricHodge())
    suite["Geometric Hodge"]["Old Form-1"] = @benchmarkable hodge_star(1, $earth, GeometricHodge())
end

begin
    suite["Inverse Diagonal Hodge"] = BenchmarkGroup()

    suite["Inverse Diagonal Hodge"]["New Form-0"] = @benchmarkable dec_inv_hodge(0, $earth, DiagonalHodge())
    suite["Inverse Diagonal Hodge"]["Old Form-0"] = @benchmarkable inv_hodge_star(0, $earth, DiagonalHodge())

    suite["Inverse Diagonal Hodge"]["New Form-1"] = @benchmarkable dec_inv_hodge(1, $earth, DiagonalHodge())
    suite["Inverse Diagonal Hodge"]["Old Form-1"] = @benchmarkable inv_hodge_star(1, $earth, DiagonalHodge())

    suite["Inverse Diagonal Hodge"]["New Form-2"] = @benchmarkable dec_inv_hodge(2, $earth, DiagonalHodge())
    suite["Inverse Diagonal Hodge"]["Old Form-2"] = @benchmarkable inv_hodge_star(2, $earth, DiagonalHodge())
end

begin
    Random.seed!(7331)
    V_1 = rand(nv(earth))
    E_1, E_2 = rand(ne(earth)), rand(ne(earth))
    T_2 = rand(ntriangles(earth))

    suite["Wedge Product"] = BenchmarkGroup()

    suite["Wedge Product"]["New Form-0, Form-1"] = @benchmarkable dec_wedge_product(Tuple{0,1}, $earth)($V_1, $E_1)
    suite["Wedge Product"]["Old Form-0, Form-1"] = @benchmarkable wedge_product(Tuple{0,1}, $earth, $V_1, $E_1)

    suite["Wedge Product"]["New Form-1, Form-1"] = @benchmarkable dec_wedge_product(Tuple{1,1}, $earth)($E_1, $E_2)
    suite["Wedge Product"]["Old Form-1, Form-1"] = @benchmarkable wedge_product(Tuple{1,1}, $earth, $E_1, $E_2)

    suite["Wedge Product"]["New Form-0, Form-2"] = @benchmarkable dec_wedge_product(Tuple{0,2}, $earth)($V_1, $T_2)
    suite["Wedge Product"]["Old Form-0, Form-2"] = @benchmarkable wedge_product(Tuple{0,2}, $earth, $V_1, $T_2)
end


# tune!(suite)

results = run(suite, verbose = true, seconds = 1)

for op in sort(collect(keys(results)))
    test = median(results[op])

    println("Operator: $op")
    for k in sort(collect(keys(test)))
        t = test[k].time / 1e6
        m = test[k].memory / 1e6
        println("Variant: $k, [$t ms, $m MB]")
    end
    println("----------------------------------------------------------------")
end