using Documenter
using Literate

@info "Loading Decapodes"
using Decapodes
using Catlab
using Catlab.WiringDiagrams
using AlgebraicPetri
using WGLMakie

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/Decapodes.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/Decapodes.jl/blob/main/docs"
end

# const literate_dir = joinpath(@__DIR__, "..", "examples")
# const generated_dir = joinpath(@__DIR__, "src", "examples")

# for (root, dirs, files) in walkdir(literate_dir)
#   out_dir = joinpath(generated_dir, relpath(root, literate_dir))
#   pmap(files) do file
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
    "Overview" => "overview.md",
    "CISM v2.1" => "cism.md",
    "Klausmeier" => "klausmeier.md",
    "Equations" => "equations.md",
    "ASCII Operators" => "ascii.md",
    "Misc Features" => "bc_debug.md",
    "Pipe Flow" => "poiseuille.md",
    "Glacial Flow" => "ice_dynamics.md",
    "Grigoriev Ice Cap" => "grigoriev.md",
    "Budyko-Sellers-Halfar" => "budyko_sellers_halfar.md",
#    "Examples" => Any[
#      "examples/cfd_example.md"
#    ],
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
