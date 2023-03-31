##############
# References #
##############

# Kotovsky, D. A. (2016), Response of the nighttime upper mesosphere to electric field changes produced by lightning discharges, Ph.D. dissertation, University of Florida, Gainesville, Florida. 

# V. P. Pasko, U. S. Inan, T. F. Bell, and Y. N. Taranenko, “Sprites produced by quasi-electrostatic heating and ionization in the lower ionosphere,” Journal of Geophysical Research: Space Physics, vol. 102, no. A3, pp. 4529–4561, 1997, doi: 10.1029/96JA03528.

# G. Veronis, V. P. Pasko, and U. S. Inan, “Characteristics of mesospheric optical emissions produced by lightning discharges,” Journal of Geophysical Research: Space Physics, vol. 104, no. A6, pp. 12645–12656, 1999, doi: 10.1029/1999JA900129.

################
# Dependencies #
################
# AlgebraicJulia
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes

# External Dependencies
using Base.MathConstants: e, π
using Distributions
using GLMakie
using Interpolations
using LinearAlgebra
using Logging
using MLStyle
using MultiScaleArrays
using OrdinaryDiffEq
using Roots
using SparseArrays

using GeometryBasics: Point2, Point3
Point2D = Point2{Float64}
Point3D = Point3{Float64}

#################
# Load the Mesh #
#################
# We assume cylindrical symmetry.
MAX_r = 400 # km
MAX_Z = 100 # km

s = loadmesh(Rectangle_30x10())
scaling_mat_to_unit_square = Diagonal([
  1/maximum(p -> p[1], s[:point]),
  1/maximum(p -> p[2], s[:point]),
  1.0])
scaling_mat_to_final_dimensions = Diagonal([
  MAX_r,
  MAX_Z,
  1.0])
scaling_mats = scaling_mat_to_final_dimensions * scaling_mat_to_unit_square
s[:point] = map(x -> scaling_mats * x, s[:point])
s[:edge_orientation] = false
orient!(s)
# Visualize the mesh.
GLMakie.wireframe(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point3D}(s)
subdivide_duals!(sd, Circumcenter())
sd[:point] = map(x -> Point3{Float64}([x[1], x[3], x[2]]), sd[:point])
sd[:dual_point] = map(x -> Point3{Float64}([x[1], x[3], x[2]]), sd[:dual_point])
sd

#############
# Constants #
#############
#temp = (return_time ./ rise_time) .^ 10 # Relative time difference.
I₀ = 250.0e3      # [kA] See Kotovsky pp. 107
T1 = 50.0         # [us] See Kotovsky Table 5-1 column 1
T2 = 1000.0       # [us]
τ₁ = T1 .*1e-6    # [s]
τ₂ = T2 .*1e-6    # [s]
ε₀ = 8.854e-12    # [F m⁻¹] Permittivity of free space
μ₀ = 4*π*1e-7     # [H m⁻¹] Permeability of free space
η = sqrt(μ₀/ε₀)   # Impedance of free space
c = sqrt(1/ε₀/μ₀) # [m s⁻¹] Speed of light
v = 2/3 * c       # [m s⁻¹] Approximate speed of propagation of lightning strike through medium
qₑ = 1.602e-19    # [coul] Charge of electron
mₑ = 9.109e-31    # [kg] Mass of electron
kB = 8.617e-5     # [eV K⁻¹] Boltzmann constant
n = 10.0          # Heidler normalizing constant
# Note that these ρ₀ and z₀ are chosen for numerical stability.
# For a rectangular grid discretization, good values are 3Δρ, 3Δz.
z₀ = 1.0e3        # [km] Heidler normalizing constant
ρ₀ = 1.0e3        # [km] Heidler normalizing constant
ρ₀ = 20.0e3        # [km] Heidler normalizing constant
a = 10.0 * 1e3    # [km] Starting altitude of vertical decay for Heidler
a = 10.0          # [m] Starting altitude of vertical decay for Heidler

constants_and_parameters = (
  qₑ = qₑ,
  mₑ = mₑ,
  ε₀ = ε₀,
  μ₀ = μ₀,
  η = η,
  kB = kB,
  c = c,
  n = n,
  e = e,
  I₀ = I₀,
  v = v,
  a = a,
  z₀ = z₀,
  ρ₀ = ρ₀)

#################
# Heidler Model #
#################
# This defines J.
# The model takes in rise time, fall time, and peak current of the lightning
# strike, (i₀).
# Note that this model does not have any tangent variables.

# See Veronis et al. § 2.
# TODO: Should this distribution should reach its max at 50 microseconds?
HeidlerForη = t -> (t / τ₁)^n / (1 + ((t /τ₁)^n)) * exp(-1.0 * t / τ₁)
ddtHeidlerForη = t -> (exp(-t / τ₁) * (t/τ₁)^n * ((t*(-(t/τ₁)^n - 1)) + n*τ₂)) /
  (t * τ₂ * ((t/τ₁)^n + 1)^2)
time_for_k = find_zero(ddtHeidlerForη, (1.0e-6, 100.0e-6))
plot(0.0:1.0e-8:100.0e-6, HeidlerForη)
a = map(0.0:1.0e-5:1.5e-3) do x
  log(HeidlerForη(x)*1e6)
end
a[1] = a[2] # Set this manually to something not -Inf
plot(0.0:1.0e-5:1.5e-3, a)
plot(0.0:1.0e-8:100.0e-6, ddtHeidlerForη)
k = HeidlerForη(time_for_k) # Risetime?

# TODO: Move this inside the Heidler Decapode.
# Perhaps we may need to use a parameter for z.
# Even though we could pre-calculate these, this would be memory intensive.
# Call this for all edges.
function compute_J(ρ, z, t)
  tret = t - z / v
  n = 10
  temp = (tret/τ₁)^n   # within parenthesis of Heidlar equation |n = 10
  I = I₀/k * temp/(1 + temp) * exp(-tret/τ₂) * (tret > 0)
  Jₛ₀ = I
  if z ≤ a
    Jₛ₀ * e ^ (-1.0 * ρ^2 / ρ₀^2)
  else
    Jₛ₀ * e ^ ((-1.0 * ρ^2 / ρ₀^2) - ((z - a)^2 / z₀^2))
  end
end

J_testing = map(sd[:point]) do p
  compute_J(p[1] *1e3, p[3] *1e3, 150.0e-6)
end
extrema(J_testing)
extrema(log.(J_testing .+ 0.001))
#mesh(s, color=J_testing)
#save("J_initial.png", mesh(s, color=J_testing))
fig_mesh, ax_mesh, ob_mesh = mesh(s, color=log.(J_testing))
ax_mesh.title = "log(J) from Heidler Model"
for t ∈ range(7.0e-6, 150.0e-6; length=400)
  sleep(0.01)
  ob_mesh.color = log.(map(sd[:point]) do p
    compute_J(p[1] *1e3, p[3] *1e3, t)
  end)
end

# Note that depending on where origin/ orientation of grid, these values may be
# pointing in the direction anti-parallel to what you expect.
# i.e. You might need to flip the sign of the argument to flatten_form.
flatten_form(vfield::Function, mesh) =  ♭(mesh,
  DualVectorField(vfield.(sd[triangle_center(sd),:dual_point])))

# Note: This is written assuming the coordinates are Cartesian.
divergence_mat = inv_hodge_star(0,sd,hodge=GeometricHodge()) * dual_derivative(1, sd)
Jₛ = flatten_form(x -> [0.0, 0.0, compute_J(x[1]*1e3, x[3]*1e3, 150.0e-6)], sd)
## Plot the divergence of Jₛ. i.e. ⋆(d(Jₛ))
mesh(s, color=divergence_mat * Jₛ, colormap=:jet)
fig_mesh, ax_mesh, ob_mesh =mesh(s, colormap=:jet, color=log.(abs.(divergence_mat *
  flatten_form(x -> [0.0, 0.0, compute_J(x[1]*1e3, x[3]*1e3, 150.0e-6)], sd))))
ax_mesh.title = "(log⋅abs⋅⋆⋅d)(J) from Heidler Model"
for t ∈ range(7.0e-6, 150.0e-6; length=400)
  sleep(0.001)
  #ob_mesh.color = map(sd[:point]) do p
  ob_mesh.color = log.(abs.(divergence_mat *
    flatten_form(x -> [0.0, 0.0, compute_J(x[1]*1e3, x[3]*1e3, t)], sd)))
end

Heidler = SummationDecapode(parse_decapode(quote
  (I, J_mask, Jₛ)::Form1{X}
  (Z, flux)::Form0{X}
  (τ₁, τ₂, I₀, v, c, n, a, k, e)::Constant{X}
  (t)::Parameter{X}
  (tret)::Form0{X}

  #I == I₀ * (one / η) * (t / τ₁)^n / (1 + (t / τ₁)^n) * e ^ (negone * t / τ₂)
  tret == t .- z ./ v
  temp == (tret./τ₁).^n
  I == I₀/k * temp ./ (1 .+ temp) .* exp(-tret ./ τ₂) .* (tret .> 0)
  Jₛ₀ == I

  J == (z .≤ a) * Jₛ₀ * exp(-1.0 .* ρ .^ 2 / ρ₀ .^ 2) +
       (z .> a) * Jₛ₀ * exp((-1.0 * ρ^2 / ρ₀^2) - ((z - a)^2 / z₀^2))
end))

#############################
# Initialize E, B, σ, and θ #
#############################

E = VForm(zeros(ne(s)))
B = DualForm{1}(zeros(ne(sd)))

# Initialize reduced electric field, θ [V m m]
θ = nothing
E₀ = nothing
# Linear model of the reduced electric field
#if option == 3
  E₀ = E
  θ = 1e21/1e6 * E ./ avg₀₁(density[:gas]) # Danny Dissert p 91 theta = E/n
#end

#Tn = zeros(nv(s)) # Tn is defined in formatAtmosphere.jl.

# Initialize conductivity
sigma = zeros(nv(s))
for p in vertices(s)
  # Note: This is assuming that points in s are in Cartesian coordinates.
  if s[p, :point][2] < 60.0
    sigma[p] = 0.0
  else
    # See Kotovsky pp.92 eq 5-2b
    # See also Pasko et al.
    sigma[p] = qe * 3.656e25 * density[:e][p] / density[:gas][p] * sqrt(200.0 / Tn[p])
  end
end

###########
# Veronis #
###########
#TODO: Handle half-timestepping in compile.
# See SymplecticIntegrators
# Primal_time: J,B,Nₑ
# Dual_time: E,Nₑ,θ,σ
# Optional flag that has symbols of TVars you want to compute.
# (i.e. this gaurantees that you don't compute things you don't want.)
# This meshes well with https://docs.sciml.ai/DiffEqDocs/latest/solvers/dynamical_solve/
# This would also come into play with DEC -> plotting toolkit.
# Tonti tried to formalize such using 3D layout of diagrams.
# (~ fields on primal-space-primal-time, or dual-space-dual-time)
# (~ particles on primal-space-dual-time, or dual-space-primal-time)
# "Dual things happening on dual steps makes sense."
# Default would be to compute primals on primal time, duals on dual time.

# Always keep in mind: B_into - B_outof = E_rightwards
# i.e. Keep your subtractions consistent.

# Assumptions that allow for cylindrical "pseudo-3D":
# - Lightning can be restricted to plane.
# - Magnetic field is orthogonal to this plane.
Veronis = SummationDecapode(parse_decapode(quote
  # TODO: Double check units on these (esp. w.r.t. σ being integrated along
  # edges.)
  B::DualForm0{X}
  E::Form1{X}
  # NOTE: You might just need σ for updating E, no need to recalculate on the
  # other time step.
  σ::Form1{X} # TODO: Maybe define this as a Form0 to avoid issue of integrating
  # conductivity.
  J::Form1{X}
  (qe,c,ε₀)::Constant{X}
  dt::Constant{X}
  θ::Form1{X}
  ρ_e::Form0{X}
  ρ_gas::Form0{X}
  Tn::Form0{X}

  # See Veronis et al. Equations 1 & 2
  Ė == -(J - σ .* E)./ε₀ + (c^2).*(⋆₁⁻¹(d̃₀(B)))
  Ė == ∂ₜ(E)

  # See Veronis et al. Equations 3
  Ḃ == ⋆₂(d₁(E))
  Ḃ == ∂ₜ(B)

  # See Kotovsky pp. 91: theta = E/n
  # Note: Updating θ with E means we are using the nonlinear model. i.e. We
  # consider electron temperature variation and species density variation.
  θ == 1e21/1e6 .* E ./ avg₀₁(ρ_gas)

  # Note: There may be a way to perform masking more efficiently,  while still
  # using the "programming language of PDEs."
  # Note: You could instead use parameters to encode things that directly
  # depend on coordinates.
  Eq5_2a_mask == θ .> (0.0603 * sqrt(200 ./ Tn))
  #Eq5_2b_mask == invert_mask(Eq5_2a_mask)
  Eq5_2b_mask == θ .≤ (0.0603 * sqrt(200 ./ Tn))
  # See Kotovsky pp. 92 5-2a
  Eq5_2a == qe * ρ_e ./ ρ_gas *
    10 .^( 50.97 + 3.026 * log10( 1e-21*θ )
      + 8.4733e-2 * log10( 1e-21*θ ).^2 )
  # See Kotovsky pp. 92 5-2b
  Eq5_2b == qe * 3.656e25 * ρ_e ./ ρ_gas .* sqrt(200.0 ./ Tn)

  σ == avg₀₁((Eq5_2a_mask .*  Eq5_2a) + (Eq5_2b_mask .*  Eq5_2b))
end))
to_graphviz(Veronis)


#############
# Operators #
#############

# TODO: Move avg₀₁ into a separate file and include, since it is not particular to this model.
function avg₀₁(sd, x)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(sd)
      append!(J, [sd[e,:∂v0],sd[e,:∂v1]])
      append!(I, [e,e])
      append!(V, [0.5, 0.5])
  end
  avg_mat = sparse(I,J,V)
  avg_mat*x
end

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :⋆₀ => x -> begin
      mat = ⋆(0,sd,hodge=hodge)
      mat*x
    end
    :⋆₁ => x -> begin
      mat = ⋆(1, sd, hodge=hodge)
      mat*x
    end
    :⋆₀⁻¹ => x -> begin 
      mat = inv_hodge_star(0, sd, hodge=hodge)
      mat*x
    end
    :⋆₁⁻¹ => x -> begin
      mat = inv_hodge_star(1, sd, hodge=hodge)
      mat*x
    end
    :d₁ => x -> begin
      mat = d(1,sd)
      mat*x
    end
    :d̃₀ => x -> begin
      mat = dual_derivative(0,sd)
      mat*x
    end
    :avg₀₁ => x -> begin
      avg_mat*x
    end
    :invert_mask => x -> (!).(x)
    :.* => (x,y) -> x .* y
    :./ => (x,y) -> x ./ y

    _ => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

#####################
# Model Composition #
#####################

# TODO

###########
# Solving #
###########

# TODO

############
# Plotting #
############

# TODO
