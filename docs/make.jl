using Documenter
using Literate
using Distributed

@info "Loading Decapodes"
using Decapodes
using Catlab
using Catlab.WiringDiagrams
using AlgebraicPetri
using CairoMakie

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/Decapodes.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/Decapodes.jl/blob/main/docs"
end

# const literate_dir = joinpath(@__DIR__, "..", "examples")
# const generated_dir = joinpath(@__DIR__, "src", "examples")

# @info "Building literate files"
# for (root, dirs, files) in walkdir(literate_dir)
#   out_dir = joinpath(generated_dir, relpath(root, literate_dir))
#   # @showprogress pmap(files) do file
#   for file in files
#     f,l = splitext(file)
#     if l == ".jl" && !startswith(f, "_")
#       Literate.markdown(joinpath(root, file), out_dir;
#         config=config, documenter=true, credit=false)
#       Literate.notebook(joinpath(root, file), out_dir;
#         execute=true, documenter=true, credit=false)
#     end
#   end
# end

@info "Building Documenter.jl docs"
makedocs(
  modules   = [Decapodes],
  format    = Documenter.HTML(
    assets = ["assets/analytics.js"],
  ),
  sitename  = "Decapodes.jl",
  doctest   = false,
  checkdocs = :none,
  pages     = Any[
    "Decapodes.jl" => "index.md",
    "Vortices" => "navier_stokes/ns.md",
    "Halfar-NS" => "halmo.md",
    "Overview" => "overview.md",
    "Klausmeier" => "klausmeier.md",
    "Glacial Flow" => "ice_dynamics.md",
    "Grigoriev Ice Cap" => "grigoriev.md",
    "Budyko-Sellers-Halfar" => "budyko_sellers_halfar.md",
    "CISM v2.1" => "cism.md",
    "NHS" => "nhs.md",
    "Cahn-Hilliard" => "cahn_hilliard.md",
    "Equations" => "equations.md",
    "ASCII Operators" => "ascii.md",
    "Misc Features" => "bc_debug.md",
    "Pipe Flow" => "poiseuille.md",
    # "Examples" => Any[
    #   "examples/cfd_example.md"
    # ],
    "Canonical Models" => "canon.md",
    "Library Reference" => "api.md"
  ]
)

@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/Decapodes.jl.git",
  branch = "gh-pages",
  devbranch = "main"
)
