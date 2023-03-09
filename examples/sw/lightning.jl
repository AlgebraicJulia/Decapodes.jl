#https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/1999JA900129

# Andrew advice:
# (Ask what operations you will want to be doing to this to determine.)
# Specify typeso f physical values and justify!
# Don't break the DEC's typing. Preserve elegance!

#TODO:
# Model Bᵩ as a Dual2Form, then take ⋆ to get a 0Form, then take d to get a
# 1Form. i.e. 

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Programs
using CombinatorialSpaces
using CombinatorialSpaces.ExteriorCalculus
using Decapodes
using MultiScaleArrays
using OrdinaryDiffEq
using MLStyle
using Distributions
using LinearAlgebra
using GLMakie
using Logging
using SparseArrays
using Base.MathConstants: e, π
# using CairoMakie 

using GeometryBasics: Point3
Point3D = Point3{Float64}
Point2D = Point2{Float64}

#################
# Heidler Model #
#################
# This defines J for us.
# Takes in rise time, fall time, and peak current of the lightning strike (i nought).

# See Veronis section 2.
# Determine normalizing constant
imax = MAX_r
dr = Δr
HeidlerForη = t -> (t / T1)^n / (1 + ((t / T1)^n)) * exp(-1.0 * t / T2)
ddtHeidlerForη = t -> (exp(-1.0 * t / T2) * (t/T1)^n * ((t*(-1.0 * (t/T1)^n - 1)) + n*T2)) /
  (t * T2 * ((t/T1)^n + 1)^2)
time_for_k = find_zero(ddtHeidlerForη, (1.0, 100.0))
k = HeidlerForη(time_for_k)
#max_index = imax/c * dr/dt; # Determine the last time index to be used
#time_vector = 0:dt:(max_index-1)*dt; # Make a vector for each time instance
I_norm_vector = ((time_vector./T1).^10)./(1+(time_vector./T1).^10).*exp(-time_vector./T2); # Calculate the outputs of the Heidlar equation without the scaling constants
#I_max = max(I_norm_vector); # Find the max value of the nonscaled heidlar equation
#I_max_index = find(I_norm_vector == I_max); # Find the index at which the max value s reached
#I_max_time = dt*(I_max_index-1); # At what time instance does the max value occur
k = I_max; # Normalization factor typically notated as eta according to danny thesis Eq 5-14
I_max =k 
η = k # Not impedance of free space.
Δz = 100.0 / 20.0
Δr = 400.0 / 30.0
z₀ = 3 * Δz # TODO: We can generalize to not use this particular grid spacing.
r₀ = 3 * Δr # TODO: We can generalize to not use this particular grid spacing.
t = 0.05
zvec = collect(0.0:Δz:100.0)
tret = t .- zvec./v # Current time instance
n = 10
temp = (tret./T1).^n   # within parenthesis of Heidlar equation |n = 10
I = Io/k * temp./(1+temp) .* exp(-tret./T2) .* (tret > 0);        
Jₛ₀ = I

vmag = 5
velocity(p) = TangentBasis(CartesianPoint(p))(vmag/4, vmag/4)
# velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, -vmag/4, 0))

#begin
u = flatten_form(velocity, earth)


# Call this for all edges.
function compute_J_mask(r, z)
  a = 10.0e3
  if z < a
    Jₛ₀ * e ^ (-1.0 * r^2 / r₀^2)
  else
    Jₛ₀ * e ^ ((-1.0 * r^2 / r₀^2) - ((z - a)^2 / z₀^2))
  end
end

Heidler = SummationDecapode(parse_decapode(quote
  (I, J_mask, Jₛ)::Form1{X} # (In the z direction)
  (Z, flux)::Form0{X} # Height (i.e. distance from ground)
  (τ₁, τ₂, I₀, v, c, one, negone, n, k, e)::Constant{X}
  (t, tret)::Parameter{X}

  #I == I₀ * (one / η) * (t / τ₁)^n / (1 + (t / τ₁)^n) * e ^ (negone * t / τ₂)
  #tret == m * dt - zvec./v # Current time instance
  tret == t - zvec./v # Current time instance
  temp == (tret./T1).^10   # within parenthesis of Heidlar equation |n = 10
  I == I₀/k * temp./(1+temp) .* exp(-tret./T2) .* (tret > 0);        
  ∂ᵨ₀()

  Jₛ == I .* J_mask
  #flux == I .* J_mask
  #Jₛ == d(flux)

end))

#############
# Constants #
#############
temp = (return_time ./ rise_time) .^ 10 # Relative time difference.
peak_current = I₀
Io = 250.0  # Pulled from Danny thesis page 107 [kA]
T1 = 50.0   # From Danny thesis Table 5-1 column 1 [us]
T2 = 1000.0 # [us]
ε₀ = 8.854e-12    # Permittivity of free space (F m-1)
μ₀ = 4*π*1e-7    # Permeability of free space (H m-1)
c = sqrt(1/ε₀/μ₀) # Speed of light

constants_and_parameters = (
  one = 1.0,
  negone = -1.0,
  qₑ = 1.602e-19,    # Charge of electron (coul)
  mₑ = 9.109e-31,    # Mass of electron (kg)
  ε₀ = ε₀,    # Permittivity of free space (F m-1)
  μ₀ = μ₀,    # Permeability of free space (H m-1)
  η = sqrt(μ₀ / ε₀), # Impedance of free space
  kB = 8.617e-5,     # Boltzmann constant (eV K-1)
  c = c, # Speed of light
  n = 10.0, # Normalizing constant for Heilder.
  e = e,
  # convert inputs to proper units
  I₀ = Io .*1e3,     # Convert from input integer to kA
  τ₁ = T1 .*1e-6,    # Convert from input integer to us
  τ₂ = T2 .*1e-6,    # Convert from input integer to us
  v = 2/3 * c        # Approximate speed of prop. of lightning strike through medium.
  )

# We assume cylindrical symmetry.
MAX_Z = 100 # km
MAX_r = 400 # km
RESOLUTION_Z = 0.2 # km
RESOLUTION_r = 0.2 # km

s = loadmesh(Rectangle_30x10())
scaling_mat_to_unit_square = Diagonal([1/maximum(p -> p[1], s[:point]),
                        1/maximum(p -> p[2], s[:point]),
                        1.0])
scaling_mat_to_final_dimensions = Diagonal([MAX_r,
                                            MAX_Z,
                                            1.0])
scaling_mats = scaling_mat_to_final_dimensions * scaling_mat_to_unit_square
s[:point] = map(x -> scaling_mats * x, s[:point])
s[:edge_orientation] = false
orient!(s)
# Visualize the mesh.
GLMakie.wireframe(s)
sd = EmbeddedDeltaDualComplex2D{Bool,Float64,Point2D}(s)
subdivide_duals!(sd, Circumcenter())
sd

#############
# Operators #
#############

# x is a Form1, y is a DualForm0 or Form2
# Returns a Form1.
function cp_2_1(x,y)
  #x_cache = t2c * x
  #y_cache = e2c * y
  x_cache = e2c * x
  y_cache = t2c * y
  broadcast!(*, y_cache, x_cache, y_cache)
  cross * y_cache
end

function edge_to_support(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      de = elementary_duals(Val{1},s, e)
      ev = CombinatorialSpaces.volume(Val{1}, s, e)
      dv = sum([dual_volume(Val{1}, s, d) for d in de])
      for d in de
          dt = incident(s, d, :D_∂e0)
          append!(I, dt)
          append!(J, fill(e, length(dt)))
          append!(V, fill(1/(dv*ev), length(dt)))
      end
  end
  sparse(I,J,V)
end

function support_to_tri(s)
  vals = Dict{Tuple{Int64, Int64}, Float64}()
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for t in 1:nv(s)
      dt = elementary_duals(Val{0},s, t)
      for d in dt
          push!(I, t)
          push!(J, d)
          push!(V, 1)
      end
  end
  sparse(I,J,V)
end

diag_vols(s) = spdiagm([dual_volume(Val{2}, s, dt) for dt in 1:nparts(s, :DualTri)])

wedge_mat(::Type{Val{(1,1)}}, s) = 2.0 * (support_to_tri(s)*diag_vols(s)*edge_to_support(s))

function pd_wedge(::Type{Val{(1,1)}}, s, α, β; wedge_t = Dict((1,1)=>wedge_mat(Val{(1,1)}, s)), kw...)
  wedge_t[(1,1)] * broadcast(*, α, β)
end

#vect(s, e) = (s[s[e,:∂v1], :point] - s[s[e,:∂v0], :point]) * sign(1, s, e)
#vect(s, e::AbstractVector) = [vect(s, el) for el in e]
#t_vects(s,t) = vect(s, triangle_edges(s,t)) .* ([1,-1,1] * sign(2,s,t))
# These helper functions convert an edge to a vector, multiplying by the sign of
# a simplex where appropriate.
edge_to_vector(s, e) = (s[e, [:∂v1, :point]] - s[e, [:∂v0, :point]]) * sign(1, s, e)
edge_to_vector(s, e::AbstractVector) = [edge_to_vector(s, el) for el in e]
t_vects(s,t) = edge_to_vector(s, triangle_edges(s,t)) .* ([1,-1,1] * sign(2,s,t))
# Return a dictionary where a key is a (triangle_idx, edge_idx) a pair, and a
# value is an index for that pair (from 1 to num_triangles*3).
# The purpose of this dictionary is so we can store values on an edge with
# respect to the triangle(s) it is a part of individually.
function comp_support(sd)
  vects = []
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(vects, (t, inds[i]))
      end
  end
  v2comp = Dict{Tuple{Int64, Int64}, Int64}()
  for (i, v) in enumerate(vects)
      v2comp[v] = i
  end
  v2comp
end
v2comp = comp_support(earth)

# Return a (num_triangles*3) by (num_edges) sparse matrix.
# Row index is the index of the (triangle_idx, edge_idx) pair according to
# comp_support.
# Column index is the index of an edge.
# Values are 1/length of the corresponding edge.
# Multiplying by this matrix thus normalizes by edge length.
function edge2comp(sd, v2comp)
  I = Vector{Int64}(); J = Vector{Int64}(); V = Vector{Float64}()
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(I, v2comp[(t,inds[i])])
          push!(J, inds[i])
          #push!(V, 1 / volume(Val{1}, sd, inds[i]))
          push!(V, 1 / SimplicialSets.volume(Val{1}, sd, inds[i]))
      end
  end
  sparse(I,J,V)
end
e2c = edge2comp(earth, v2comp)

# Return a (num_triangles*3) by (num_triangles) sparse matrix.
# Row index is the index of the (triangle_idx, edge_idx) pair according to
# comp_support.
# Column index is the index of a triangle.
# Values are 1.
function tri2comp(sd, v2comp)
  I = Vector{Int64}(); J = Vector{Int64}(); V = Vector{Float64}()
  for t in 1:ntriangles(sd)
      inds = triangle_edges(sd, t)
      for i in 1:3
          push!(I, v2comp[(t,inds[i])])
          push!(J, t)
          push!(V, 1)
      end
  end
  sparse(I,J,V)
end
t2c = tri2comp(earth, v2comp)

# Return a (num_edges) by (num_triangles*3) sparse matrix.
function changes(s, v2comp)
  orient_vals = [1,-1,1]
  I = Vector{Int64}(); J = Vector{Int64}(); V = Vector{Float64}()
  for t in 1:ntriangles(s)
    inds = triangle_edges(s, t)
    e_vects = t_vects(s,t)
    for i in 1:3
      ns = [(i+1)%3 + 1, i%3+1]
      ort = e_vects[i] × (e_vects[i] × e_vects[ns[1]])
      n_ort = normalize(ort)
      append!(I, inds[ns[1]])
      append!(J, v2comp[(t,inds[i])])
      append!(V, dot(n_ort, e_vects[ns[1]]) * orient_vals[ns[1]] * sign(1, s, ns[1])* orient_vals[i]* sign(2,s,t) / 3.0)
      append!(I, inds[ns[2]])
      append!(J, v2comp[(t,inds[i])])
      append!(V, dot(n_ort, e_vects[ns[2]]) * orient_vals[ns[2]] * sign(1, s, ns[2])* orient_vals[i]* sign(2,s,t) / 3.0)
    end
  end
  sparse(I,J,V, ne(s), ntriangles(s)*3)
end
cross = changes(earth, v2comp)

function avg_mat(::Type{Val{(0,1)}},s)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(s)
      append!(J, [s[e,:∂v0],s[e,:∂v1]])
      append!(I, [e,e])
      append!(V, [0.5, 0.5])
  end
  sparse(I,J,V)
end

function avg₀₁(sd, x)
  I = Vector{Int64}()
  J = Vector{Int64}()
  V = Vector{Float64}()
  for e in 1:ne(sd)
      append!(J, [sd[e,:∂v0],sd[e,:∂v1]])
      append!(I, [e,e])
      append!(V, [0.5, 0.5])
  end
  sparse(I,J,V)*x
end

hodge = GeometricHodge()
d₀(x,sd) = d(0,sd)*x
d₁(x,sd) = d(1,sd)*x
⋆₀(x,sd,hodge) = ⋆(0,sd,hodge=hodge)*x
⋆₁(x,sd,hodge) = ⋆(1,sd,hodge=hodge)*x
⋆₂(x,sd,hodge) = ⋆(2,sd,hodge=hodge)*x
⋆₀⁻¹(x,sd,hodge) = inv_hodge_star(0,sd,hodge=hodge)*x
⋆₁⁻¹(x,sd,hodge) = inv_hodge_star(1,sd,hodge=hodge)*x
∧₁₁′(x,y,sd) = pd_wedge(Val{(1,1)}, sd, x, y)
∧₁₀′(x,y) = cp_2_1(x,y)
i₀′(x,y,sd,hodge) = -1.0 * (∧₁₀′(x, ⋆₂(y,sd,hodge)))
i₁′(x,y,sd,hodge) = ⋆₀⁻¹(∧₁₁′(x, ⋆₁(y,sd,hodge), sd),sd,hodge) #⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
L1′(x,y,sd,hodge) = i₀′(x,d₁(y,sd),sd,hodge) + d₀(i₁′(x,y,sd,hodge),sd)

function generate(sd, my_symbol; hodge=GeometricHodge())
  i0 = (v,x) -> ⋆(1, sd, hodge=hodge)*wedge_product(Tuple{0,1}, sd, v, inv_hodge_star(0,sd, hodge=DiagonalHodge())*x)
  op = @match my_symbol begin
    :k => x->2000x
    :μ => x->-0.0001x
    # :μ => x->-2000x
    :β => x->2000*x
    :γ => x->1*x
    :⋆₀ => x->⋆(0,sd,hodge=hodge)*x
    :⋆₁ => x->⋆(1, sd, hodge=hodge)*x
    :⋆₀⁻¹ => x->inv_hodge_star(0,sd, x; hodge=hodge)
    :⋆₁⁻¹ => x->inv_hodge_star(1,sd,hodge=hodge)*x
    :d₀ => x->d(0,sd)*x
    :d₁ => x->d(1,sd)*x
    :dual_d₀ => x->dual_derivative(0,sd)*x
    :dual_d₁ => x->dual_derivative(1,sd)*x
    :d̃₁ => x->dual_derivative(1,sd)*x
    :∧₀₁ => (x,y)-> wedge_product(Tuple{0,1}, sd, x, y)
    :∧₁₀ => (x,y)-> wedge_product(Tuple{1,0}, sd, x, y)
    :plus => (+)
    :(-) => x-> -x
    # :L₀ => (v,x)->dual_derivative(1, sd)*(i0(v, x))
    :L₀ => (v,x)->begin
      dual_derivative(1, sd)*⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x) # wedge product of 1 and 2 form???  switched with 183; add in lv*d term
      # ⋆(1, sd; hodge=hodge)*wedge_product(Tuple{1,0}, sd, v, x)
    end
    :i₀ => i0 
    :Δ₀ => x -> begin # dδ
      δ(1, sd, d(0, sd)*x, hodge=hodge) end
    # :Δ₀ => x -> begin # d ⋆ d̃ ⋆⁻¹
    #   y = dual_derivative(1,sd)*⋆(1, sd, hodge=hodge)*d(0,sd)*x
    #   inv_hodge_star(0,sd, y; hodge=hodge)
    # end
    :Δ₁ => x -> begin # dδ + δd
      δ(2, sd, d(1, sd)*x, hodge=hodge) + d(0, sd)*δ(1, sd, x, hodge=hodge)
    end

    :δ₁ => x -> inv_hodge_star(0, sd, hodge=hodge) * dual_derivative(1,sd) * ⋆(1, sd, hodge=hodge) * x
    #:i₁′ => (v,x) -> inv_hodge_star(0,sd, hodge=hodge) * wedge_product(Tuple{1,1}, sd, v, ⋆(1, sd, hodge=hodge) * x) #⋆₀⁻¹{X}(∧₁₁′(F1, ⋆₁{X}(F1′)))
    :i₁′ => (x,y) -> i₁′(x,y,sd,hodge)
    #:L₁′ = ... + d(0,sd)*i₁′(v,x) #i₀′(F1, d₁{X}(F1′)) + d₀{X}(i₁′(F1, F1′))
    :L₁′ => (x,y) -> L1′(x,y,sd,hodge)
    :neg₁ => x -> -1.0 * x
    :neg₀ => x -> -1.0 * x
    :half => x -> 0.5 * x
    :third => x -> x / 3.0
    #:R₀ => x-> 1.38064852e-23 * 6.0221409e23 / (28.96 / 1000) # Boltzmann constant * ??? / (molecular mass / 1000)
    :R₀ => x-> boltzmann_constant * 6.0221409e23 / (mol_mass / 1000) # Boltzmann constant * ??? / (molecular mass / 1000)
    # :kᵥ => x->0.000210322*x # / density
    :kᵥ => x->1.2e-5*x # / density
    # These are the steps used to compute k.
    # We have no boundaries, so I set k to the constant k₁
    #kₜ = 0.0246295028571 # Thermal conductivity
    #k_cyl = kₜ * 4
    #density = 0.000210322
    #cₚ = 1004.703 # Specific Heat at constant pressure
    #k₁ = kₜ / (density * cₚ) # Heat diffusion constant in fluid
    #k₂ = k_cyl / (density * cₚ) # Heat diffusion constant in cylinder
    #k_col = fill(k₁, ne(s))
    #k_col[cyl] .= k₂
    #k = diagm(k_col)
    :k => x->k₁*x
    :div₀ => (v,x) -> v / x
    :div₁ => (v,x) -> v / x
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(sd)
          append!(J, [sd[e,:∂v0],sd[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      sparse(I,J,V)*x
    end
    :.* => (x,y) -> x .* y
    :./ => (x,y) -> x ./ y

    # :Δ₁ => x -> begin # d ⋆ d̃ ⋆⁻¹ + ⋆ d̃ ⋆ d
    #   y = dual_derivative(0,sd)*⋆(2, sd, hodge=hodge)*d(1,sd)*x
    #   inv_hodge_star(2,sd, y; hodge=hodge) 
    #   z = d(0, sd) * inverse_hode_star(2, sd, dual_derivative(1, sd)*⋆(1,sd, hodge=hodge)*x; hodge=hodge)
    #   return y + z
    # end
    :debug => (args...)->begin println(args[1], length.(args[2:end])) end
    x=> error("Unmatched operator $my_symbol")
  end
  # return (args...) -> begin println("applying $my_symbol"); println("arg length $(length.(args))"); op(args...);end
  return (args...) ->  op(args...)
end

flatten_form(vfield::Function, mesh) =  ♭(mesh, DualVectorField(vfield.(mesh[triangle_center(mesh),:dual_point])))

##########################
# Constants & Parameters #
##########################


######################
# Initial Conditions #
######################

begin
  #vmag = 500
  vmag = 5
  # velocity(p) = vmag*ϕhat(p)
  velocity(p) = TangentBasis(CartesianPoint(p))(vmag/4, vmag/4)
  # velocity(p) = TangentBasis(CartesianPoint(p))((vmag/4, -vmag/4, 0))

# visualize the vector field
  ps = earth[:point]
  #ns = ((x->x) ∘ (x->Vec3f(x...))∘velocity).(ps)
  #GLMakie.arrows(
  #    ps, ns, fxaa=true, # turn on anti-aliasing
  #    linecolor = :gray, arrowcolor = :gray,
  #    linewidth = 20.1, arrowsize = 20*Vec3f(3, 3, 4),
  #    align = :center, axis=(type=Axis3,)
  #)
end

#begin
u = flatten_form(velocity, earth)
#c_dist = MvNormal([RADIUS/√(2), RADIUS/√(2)], 20*[1, 1])
#c = 100*[pdf(c_dist, [p[1], p[2]]) for p in earth[:point]]

theta_start = 45*pi/180
phi_start = 0*pi/180
x = RADIUS*cos(phi_start)*sin(theta_start)
y = RADIUS*sin(phi_start)*sin(theta_start)
z = RADIUS*cos(theta_start)
#c_dist₁ = MvNormal([x, y, z], 20*[1, 1, 1])
#c_dist₂ = MvNormal([x, y, -z], 20*[1, 1, 1])
#c_dist = MixtureModel([c_dist₁, c_dist₂], [0.6,0.4])

#pfield = 100000*[abs(p[3]) for p in earth[:point]]

# TODO What are good initial conditions for this?
#n = 100*[pdf(c_dist, [p[1], p[2], p[3]]) for p in earth[:point]]
#n = 100000*[p[3] for p in earth[:point]]
# Observe: Set n of protons out of phase with n of N2.
N2_n = 100_000*[p[3] > 0.0 ? abs(p[3]) : 0.0 for p in earth[:point]]
proton_n = 100_000*[p[3] < 0.0 ? abs(p[3]) : 0.0 for p in earth[:point]]

temperature = 300 # kelvin
N2_pfield = N2_n * R₀ * temperature
proton_pfield = proton_n * R₀ * temperature

#m = 8e22 # dipole moment units: A*m^2
#μ₀ = 4*π*1e-7 # permeability units: kg*m/s^2/A^2
## Bθ(p) = -μ₀*m*sin(theta(p))/(4*π*RADIUS^3)  # RADIUS instead of r(p)*1000 
## B1Form = flatten_form(p->TangentBasis(CartesianPoint(p))(Bθ(CartesianPoint(p)), 0), earth)
#Br(p) = -μ₀*2*m*cos(theta(p))/(4*π*RADIUS^3)  # RADIUS instead of r(p)*1000 
#Br_flux = hodge_star(earth, TriForm(map(triangles(earth)) do t 
#                        dual_pid = triangle_center(earth, t)
#                        p = earth[dual_pid, :dual_point]
#                        return Br(CartesianPoint(p))
#                        end))
## B₀(p) = sqrt(Bθ(p)^2+Br(p)^2)

###########
# Solving #
###########

#@info("Precompiling Solver")
#fₘ(Nothing, u₀, my_constants, (0, 1e-8))
#prob = ODEProblem(fₘ,u₀,(0,1e-4),my_constants)
#soln = solve(prob, Tsit5(), progress=true)
#soln.retcode != :Unstable || error("Solver was not stable")
#@info("Solving")
#prob = ODEProblem(fₘ,u₀,(0,tₑ),my_constants)
#soln = solve(prob, Tsit5())
#@info("Done")
#end

#begin
#mass(soln, t, mesh, concentration=:P) = sum(⋆(0, mesh)*findnode(soln(t), concentration))
#
#@show extrema(mass(soln, t, earth, :P) for t in 0:tₑ/150:tₑ)
#end
#mesh(primal_earth, color=findnode(soln(0), :P), colormap=:jet)
#mesh(primal_earth, color=findnode(soln(0) - soln(tₑ), :P), colormap=:jet)
#begin

############
# Plotting #
############

# Plot the result
times = range(0.0, tₑ, length=150)
colors_proton = [findnode(soln(t), :proton_n) for t in times]
colors_N2 = [findnode(soln(t), :N2_n) for t in times]
# Initial frame
fig = GLMakie.Figure()
p1 = GLMakie.mesh(fig[1,2], primal_earth, color=colors_proton[1], colormap=:jet, colorrange=extrema(colors_proton[1]))
p2 = GLMakie.mesh(fig[1,3], primal_earth, color=colors_N2[1], colormap=:jet, colorrange=extrema(colors_N2[1]))
Colorbar(fig[1,1], ob_proton)
Colorbar(fig[1,4], ob_N2)
Label(fig[1,2,Top()], "Proton n")
Label(fig[1,3,Top()], "N2 n")
lab1 = Label(fig[1,2,Bottom()], "")
lab2 = Label(fig[1,3,Bottom()], "")

# Animation
using Printf
record(fig, "lightning.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    p1.plot.color = findnode(soln(t), :proton_n)
    p2.plot.color = findnode(soln(t), :N2_n)
    lab1.text = @sprintf("%.2f",t)
    lab2.text = @sprintf("%.2f",t)
end

times = range(0.0, tₑ, length=150)
colors_proton = [findnode(soln(t), :proton_P) for t in times]
colors_N2 = [findnode(soln(t), :N2_P) for t in times]
# Initial frame
fig = GLMakie.Figure()
p1 = GLMakie.mesh(fig[1,2], primal_earth, color=colors_proton[1], colormap=:jet, colorrange=extrema(colors_proton[1]))
p2 = GLMakie.mesh(fig[1,3], primal_earth, color=colors_N2[1], colormap=:jet, colorrange=extrema(colors_N2[1]))
Colorbar(fig[1,1], ob_proton)
Colorbar(fig[1,4], ob_N2)
Label(fig[1,2,Top()], "Proton P")
Label(fig[1,3,Top()], "N2 P")
lab1 = Label(fig[1,2,Bottom()], "")
lab2 = Label(fig[1,3,Bottom()], "")

# Animation
using Printf
record(fig, "lightning.gif", range(0.0, tₑ; length=150); framerate = 30) do t
    p1.plot.color = findnode(soln(t), :proton_P)
    p2.plot.color = findnode(soln(t), :N2_P)
    lab1.text = @sprintf("%.2f",t)
    lab2.text = @sprintf("%.2f",t)
end
#end

########################
# Interactive Plotting #
########################

#function interactive_sim_view(my_mesh::EmbeddedDeltaSet2D, tₑ, soln; loop_times = 1)
#  times = range(0.0, tₑ, length = 150)
#  colors = [findnode(soln(t), :n) for t in times]
#  fig, ax, ob = GLMakie.mesh(my_mesh, color=colors[1],
#    colorrange = extrema(colors[1]), colormap=:jet)
#  display(fig)
#  loop = range(0.0, tₑ; length=150)
#  for _ in 1:loop_times
#    for t in loop
#      ob.color = findnode(soln(t), :n)
#      sleep(0.05)
#    end
#    for t in reverse(loop)
#      ob.color = findnode(soln(t), :n)
#      sleep(0.05)
#    end
#  end
#end
#
#interactive_sim_view(primal_earth, tₑ, soln, loop_times = 10)
