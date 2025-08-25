using CombinatorialSpaces

s = EmbeddedDeltaSet3D{Bool,Point3d}()
add_vertices!(s, 4, point=[Point3d(0,0,0), Point3d(1,0,0),
  Point3d(0,1,0), Point3d(0,0,1)])
glue_tetrahedron!(s, 1, 2, 3, 4)
orient!(s)
sd = EmbeddedDeltaDualComplex3D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

wdg_01 = dec_wedge_product(Tuple{0,1}, sd)
star₁_mat = hodge_star(1,sd,DiagonalHodge())
dual_d₂_mat = dual_derivative(2,sd)
inv_star₃_mat = inv_hodge_star(0,sd,DiagonalHodge())

codif = inv_star₃_mat * dual_d₂_mat * star₁_mat

T = Float64[1,0,0,0]
u = eval_constant_primal_form(sd, Point3d(1,0,0))

@show -(codif * wdg_01(T, u))

s = EmbeddedDeltaSet2D{Bool,Point3d}()
add_vertices!(s, 3, point=[Point3d(0,0,0), Point3d(1,0,0),
  Point3d(0,1,0)])
glue_triangle!(s, 1, 2, 3)
orient!(s)
sd = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s)
subdivide_duals!(sd, Circumcenter())

wdg_01 = dec_wedge_product(Tuple{0,1}, sd)
star₁_mat = hodge_star(1,sd,DiagonalHodge())
dual_d₁_mat = dual_derivative(1,sd)
inv_star₃_mat = inv_hodge_star(0,sd,DiagonalHodge())

codif = inv_star₃_mat * dual_d₁_mat * star₁_mat

T = Float64[1,0,0]
u = eval_constant_primal_form(sd, Point3d(1,0,0))

@show -(codif * wdg_01(T, u))