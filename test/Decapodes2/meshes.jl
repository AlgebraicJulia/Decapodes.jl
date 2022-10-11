using Catlab, Catlab.CategoricalAlgebra
using CombinatorialSpaces
using Test

# TODO: Remove this line when we can do using Decapodes.
include("../../src/Decapodes2/meshes.jl")

# TODO: Move the testset macro to runtests.jl.
@testset "Meshes" begin

# Unit Icosphere
#unit_icosphere = loadmesh(UnitIcosphere())
#@test nv(unit_icosphere) == 2562
#@test ne(unit_icosphere) == 7680
#@test nparts(unit_icosphere, :Tri) == 5120
#
## Icosphere at minimum altitude of thermosphere
#thermo_icosphere = loadmesh(ThermoIcosphere())
#@test nv(thermo_icosphere) == 2562
#@test ne(thermo_icosphere) == 7680
#@test nparts(thermo_icosphere, :Tri) == 5120
#
## Unit UV sphere (TIE GCM grid)
#unit_uvsphere = loadmesh(UnitUVSphere())
#@test nv(unit_uvsphere) == 2522
#@test ne(unit_uvsphere) == 7560
#@test nparts(unit_uvsphere, :Tri) == 5040
#
## Icosphere at minimum altitude of thermosphere
#thermo_uvsphere = loadmesh(ThermoUVSphere())
#@test nv(thermo_uvsphere) == 2522
#@test ne(thermo_uvsphere) == 7560
#@test nparts(thermo_uvsphere, :Tri) == 5040
#
#
#magnitude = (sqrt ∘ (x -> foldl(+, x*x)))
#
## The definition of a discretization of a sphere of unspecified radius.
#ρ′ = magnitude(unit_icosphere[:point][begin])
#@test all(isapprox.(magnitude.(unit_icosphere[:point]), ρ′))
## The definition of a discretization of a sphere of radius ρ.
#ρ = 1
#@test all(isapprox.(magnitude.(unit_icosphere[:point]), ρ))
#
## The definition of a discretization of a sphere of unspecified radius.
#ρ′ = magnitude(thermo_icosphere[:point][begin])
#@test all(isapprox.(magnitude.(thermo_icosphere[:point]), ρ′))
#ρ = 6371+90
## The definition of a discretization of a sphere of radius ρ.
#@test all(isapprox.(magnitude.(thermo_icosphere[:point]), ρ))
#
## The definition of a discretization of a sphere of unspecified radius.
#ρ′ = magnitude(unit_uvsphere[:point][begin])
#@test all(isapprox.(magnitude.(unit_uvsphere[:point]), ρ′))
## The definition of a discretization of a sphere of radius ρ.
#ρ = 1
#@test all(isapprox.(magnitude.(unit_uvsphere[:point]), ρ))
#
## The definition of a discretization of a sphere of unspecified radius.
#ρ′ = magnitude(thermo_uvsphere[:point][begin])
#@test all(isapprox.(magnitude.(thermo_uvsphere[:point]), ρ′))
## The definition of a discretization of a sphere of radius ρ.
#ρ = 6371+90
#@test all(isapprox.(magnitude.(thermo_uvsphere[:point]), ρ))

end
