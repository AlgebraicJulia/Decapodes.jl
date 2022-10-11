using Catlab, Catlab.CategoricalAlgebra
using CombinatorialSpaces
using Test

# TODO: Remove this line when we can do using Decapodes.
include("../../src/Decapodes2/meshes.jl")

# TODO: Move the testset macro to runtests.jl.
@testset "Meshes" begin

magnitude = (sqrt ∘ (x -> foldl(+, x*x)))
unit_radius = 1

# Unit Icospheres are of the proper dimensions, and are spheres.
unit_icosphere1 = loadmesh(UnitIcosphere1())
@test nv(unit_icosphere1) == 12
@test ne(unit_icosphere1) == 30
@test nparts(unit_icosphere1, :Tri) == 20
ρ = magnitude(unit_icosphere1[:point][begin])
@test all(isapprox.(magnitude.(unit_icosphere1[:point]), ρ))
@test ρ == unit_radius

unit_icosphere2 = loadmesh(UnitIcosphere2())
@test nv(unit_icosphere2) == 42
@test ne(unit_icosphere2) == 120
@test nparts(unit_icosphere2, :Tri) == 80
ρ = magnitude(unit_icosphere2[:point][begin])
@test all(isapprox.(magnitude.(unit_icosphere2[:point]), ρ))
@test ρ == unit_radius

unit_icosphere3 = loadmesh(UnitIcosphere3())
@test nv(unit_icosphere3) == 162
@test ne(unit_icosphere3) == 480
@test nparts(unit_icosphere3, :Tri) == 320
ρ = magnitude(unit_icosphere3[:point][begin])
@test all(isapprox.(magnitude.(unit_icosphere3[:point]), ρ))
@test ρ == unit_radius

unit_icosphere4 = loadmesh(UnitIcosphere4())
@test nv(unit_icosphere4) == 642
@test ne(unit_icosphere4) == 1920
@test nparts(unit_icosphere4, :Tri) == 1280
ρ = magnitude(unit_icosphere4[:point][begin])
@test all(isapprox.(magnitude.(unit_icosphere4[:point]), ρ))
@test ρ == unit_radius

unit_icosphere5 = loadmesh(UnitIcosphere5())
@test nv(unit_icosphere5) == 2562
@test ne(unit_icosphere5) == 7680
@test nparts(unit_icosphere5, :Tri) == 5120
ρ = magnitude(unit_icosphere5[:point][begin])
@test all(isapprox.(magnitude.(unit_icosphere5[:point]), ρ))
@test ρ == unit_radius

end

