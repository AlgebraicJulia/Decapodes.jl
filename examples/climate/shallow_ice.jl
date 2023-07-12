# AlgebraicJulia Dependencies
using Catlab
using Catlab.Graphics
using CombinatorialSpaces
using Decapodes

# External Dependencies
using MLStyle
using MultiScaleArrays
using LinearAlgebra
using OrdinaryDiffEq
using JLD2
using SparseArrays
using GLMakie
using GeometryBasics: Point2
Point2D = Point2{Float64}

####################
# Define the model #
####################

# Equation 2 from Halfar, P. ON THE DYNAMICS OF THE ICE SHEETS. (1981)
halfar_eq2 = @decapode begin
  (h,Γ)::Form0
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ .* d(h) .* abs(d(h))^(n-1) .* avg₀₁(h^(n+2)))
end
to_graphviz(halfar_eq2)

# Equation 1 from Glen, J. W. THE FLOW LAW OF ICE: A discussion of the
# assumptions made in glacier theory, their experimental foundations and
# consequences. (1958)
glens_law = @decapode begin
  Γ::Form0
  (A,ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end
to_graphviz(glens_law)

ice_dynamics_composition_diagram = @relation () begin
  dynamics(h,Γ,n)
  stress(Γ,n)
end
to_graphviz(ice_dynamics_composition_diagram, box_labels=:name, junction_labels=:variable, prog="circo")

ice_dynamics_cospan = oapply(ice_dynamics_composition_diagram,
  [Open(halfar_eq2, [:h,:Γ,:n]),
  Open(glens_law, [:Γ,:n])])

ice_dynamics = apex(ice_dynamics_cospan)
to_graphviz(ice_dynamics)

ice_dynamics = expand_operators(ice_dynamics)
to_graphviz(ice_dynamics)

infer_types!(ice_dynamics)
to_graphviz(ice_dynamics)

resolve_overloads!(ice_dynamics)
to_graphviz(ice_dynamics)

# Demonstrate storing as JSON.
write_json_acset(ice_dynamics, "ice_dynamics.json")
# When reading back in, we specify that all attributes are "Symbol"s.
ice_dynamics2 = read_json_acset(SummationDecapode{Symbol,Symbol,Symbol}, "ice_dynamics.json")
to_graphviz(ice_dynamics2)
# Or, you could choose to interpret the data as "String"s.
ice_dynamics3 = read_json_acset(SummationDecapode{String,String,String}, "ice_dynamics.json")
to_graphviz(ice_dynamics3)

###################
# Define the mesh #
###################

include("../../grid_meshes.jl")
s′ = triangulated_grid(10_000,10_000,800,800,Point3D)
s = EmbeddedDeltaDualComplex2D{Bool, Float64, Point3D}(s′)
subdivide_duals!(s, Barycenter())

########################################################
# Define constants, parameters, and initial conditions #
########################################################

n = 3
ρ = 910
g = 9.8
A = 1e-16
Γ = (2/(n+2))*A*(ρ*g)^n

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s′)) do (x,y)
  10 + (1e-7)*((x-32)^2 + (y-32)^2)
end

# Visualize initial condition for ice sheet height.
mesh(s′, color=h₀, colormap=:jet)

# Store these values to be passed to the solver.
u₀ = construct(PhysicsState, [VectorForm(h₀)], Float64[], [:h])
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g,
  stress_A = A)

#############################################
# Define how symbols map to Julia functions #
#############################################

hodge = GeometricHodge()
function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :d₀ => x -> begin
      d₀ = d(s,0)
      d₀ * x
    end
    :dual_d₀ => x -> begin
      dual_d₀ = dual_derivative(s,0)
      dual_d₀ * x
    end
    :⋆₁ => x -> begin
      ⋆₁ = ⋆(s,1)
      ⋆₁ * x
    end
    :⋆₀⁻¹ => x -> begin
      ⋆₀⁻¹ = inv_hodge_star(s,0)
      ⋆₀⁻¹ * x
    end
    :avg₀₁ => x -> begin
      I = Vector{Int64}()
      J = Vector{Int64}()
      V = Vector{Float64}()
      for e in 1:ne(s)
          append!(J, [s[e,:∂v0],s[e,:∂v1]])
          append!(I, [e,e])
          append!(V, [0.5, 0.5])
      end
      avg_mat = sparse(I,J,V)
      avg_mat * x
    end
    :^ => (x,y) -> x .^ y
    :* => (x,y) -> x .* y
    :abs => x -> abs.(x)
    :show => x -> begin
      println(x)
      x
    end
    x => error("Unmatched operator $my_symbol")
  end
  return (args...) -> op(args...)
end

#######################
# Generate simulation #
#######################

sim = eval(gensim(ice_dynamics, dimension=2))
fₘ = sim(s, generate)

##################
# Run simulation #
##################

tₑ = 3e5

# Julia will pre-compile the generated simulation the first time it is run.
@info("Precompiling Solver")
prob = ODEProblem(fₘ, u₀, (0, 1e-8), constants_and_parameters)
soln = solve(prob, Tsit5())
soln.retcode != :Unstable || error("Solver was not stable")

# This next run should be fast.
@info("Solving")
prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
soln = solve(prob, Tsit5())
@show soln.retcode
@info("Done")

@save "ice_dynamics.jld2" soln

#############
# Visualize #
#############

# Visualize the final conditions.
mesh(s′, color=findnode(soln(tₑ), :h), colormap=:jet)

# Create a gif
begin
  frames = 100
  fig, ax, ob = GLMakie.mesh(s′, color=findnode(soln(0), :h), colormap=:jet, colorrange=extrema(findnode(soln(tₑ), :h)))
  Colorbar(fig[1,2], ob)
  record(fig, "ice_dynamics.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
    ob.color = findnode(soln(t), :h)
  end
end


begin end









## Equation 1 from Halfar, P. ON THE DYNAMICS OF THE ICE SHEETS. (1981)
#halfar_eq1 = @decapode begin
#  h::Form0
#  (Γ,n)::Constant
#
#  ḣ == ∂ₜ(h)
#  #ḣ == -1 * Γ * ∘(⋆, d, ⋆)(d(h) .* abs(d(h))^(n-1) .* avg₀₁(h^(n+2)))
#  ḣ == Γ * ∘(⋆, d, ⋆)(d(h) .* abs(d(h))^(n-1) .* avg₀₁(h^(n+2)))
#  #ḣ == ∘(⋆, d, ⋆)(Γ * avg₀₁(h^(n+2)) .* abs(d(h))^(n-1) .* d(h))
#end
#
## Equation 1 from Glen, J. W. THE FLOW LAW OF ICE: A discussion of the
## assumptions made in glacier theory, their experimental foundations and
## consequences. (1958)
#glen_eq1 = @decapode begin
#  (ϵ,σ,A,n)::Constant
#
#  ϵ == (σ/A)^n
#end
#
## Equation 1 from Mahaffy, M. W. A THREE-DIMENSIONAL NUMERICAL MODEL OF ICE
## SHEETS: Tests on the Barnes ICe Cap, Northwest Territories. (1976)
#mahaffy_eq1 = @decapode begin
#  h::Form0
#  (ϵ,σ,A,n)::Constant
#  b::Parameter
#
#  ḣ == ∂ₜ(h)
#  ḣ == b - ∘(⋆, d, ⋆)(q)
#end
#
## Equations 9 and 10 from Mahaffy, M. W. A THREE-DIMENSIONAL NUMERICAL MODEL OF
## ICE SHEETS: Tests on the Barnes ICe Cap, Northwest Territories. (1976)
#mahaffy_eqs9_10 = @decapode begin
#  h::Form0
#  (σ,ρ,g)::Constant
#
#  σ == -ρ * g * (h-z) * d(h)
#end
#
## Equation 13 from Mahaffy, M. W. A THREE-DIMENSIONAL NUMERICAL MODEL OF ICE
## SHEETS: Tests on the Barnes ICe Cap, Northwest Territories. (1976)
#mahaffy_eq13 = @decapode begin
#  h::Form0
#  α::Form1
#
#  α == abs(d(h))
#end
#
## Equations 14 and 15 from Mahaffy, M. W. A THREE-DIMENSIONAL NUMERICAL MODEL OF
## ICE SHEETS: Tests on the Barnes ICe Cap, Northwest Territories. (1976)
#mahaffy_eqs14_15 = @decapode begin
#  h::Form0
#  (α,q)::Form1
#
#end
## Infer the forms of dependent variables, and resolve which versions of DEC
## operators to use.
#halfar = expand_operators(halfar_eq1)
#infer_types!(expand_operators(ice_dynamics))
#resolve_overloads!(ice_dynamics)
#infer_types!(halfar, op1_inf_rules_2D,
#  vcat(op2_inf_rules_2D, [
#    (proj1_type = :Constant, proj2_type = :Literal, res_type = :Constant, op_names = [:/, :./, :*, :.*]),
#    (proj1_type = :Literal, proj2_type = :Constant, res_type = :Constant, op_names = [:/, :./, :*, :.*])]
#  ))
#resolve_overloads!(halfar)
#
## Visualize.
#to_graphviz(halfar)
#to_graphviz(ice_dynamics)
#constants_and_parameters = (
#  n = n,
#  ρ = ρ,
#  g = g,
#  A = A,
#  Γ = Γ)
#prob = ODEProblem(fₘ, u₀, (0, tₑ), constants_and_parameters)
#soln = solve(prob, Tsit5())
#extrema(findnode(soln(0.0), :h))
#extrema(findnode(soln(tₑ ), :h))
#
#extrema(d0s * findnode(soln(0.0), :h))
#extrema(d0s * findnode(soln(tₑ ), :h))

#sim = eval(gensim(halfar, dimension=2))
#fₘ = sim(s, generate)
#@save "halfar.jld2" soln
#begin # BEGIN Gif creation
#frames = 100
## Initial frame
##fig = GLMakie.Figure(resolution = (400, 400))
#fig, ax, ob = GLMakie.mesh(s′, color=findnode(soln(0), :h), colormap=:jet, colorrange=extrema(findnode(soln(tₑ), :h)))
#Colorbar(fig[1,2], ob)
#
## Animation
#record(fig, "halfar.gif", range(0.0, tₑ; length=frames); framerate = 30) do t
#    ob.color = findnode(soln(t), :h)
#end
#
#end # END Gif creation
