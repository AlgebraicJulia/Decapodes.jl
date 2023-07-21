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
  h::Form0
  Γ::Form1
  n::Constant

  ḣ == ∂ₜ(h)
  ḣ == ∘(⋆, d, ⋆)(Γ * d(h) * avg₀₁(mag(♯(d(h)))^(n-1)) * avg₀₁(h^(n+2)))
end
to_graphviz(halfar_eq2)

# Equation 1 from Glen, J. W. THE FLOW LAW OF ICE: A discussion of the
# assumptions made in glacier theory, their experimental foundations and
# consequences. (1958)
glens_law = @decapode begin
  (A, Γ)::Form1
  (ρ,g,n)::Constant
  
  Γ == (2/(n+2))*A*(ρ*g)^n
end
to_graphviz(glens_law)

#####################
# Compose the model #
#####################

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

###################
# Store the model #
###################

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
A = fill(1e-16, ne(s))

# Ice height is a primal 0-form, with values at vertices.
h₀ = map(point(s′)) do (x,y)
  (7072-((x-5000)^2 + (y-5000)^2)^(1/2))/9e3+10
end

# Visualize initial condition for ice sheet height.
mesh(s′, color=h₀, colormap=:jet)

# Store these values to be passed to the solver.
u₀ = construct(PhysicsState, [VectorForm(h₀), VectorForm(A)], Float64[], [:h, :stress_A])
constants_and_parameters = (
  n = n,
  stress_ρ = ρ,
  stress_g = g)

#############################################
# Define how symbols map to Julia functions #
#############################################

function generate(sd, my_symbol; hodge=GeometricHodge())
  op = @match my_symbol begin
    :♯ => x -> begin
      ♯(sd, EForm(x))
  end
    :mag => x -> begin
      norm.(x)
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

tₑ = 5e13

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
