@info "Loading Documenter"
using Documenter
using Literate
using Distributed

@info "Loading Decapodes"
using Decapodes

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
  draft = true;
  pagesonly = true,
  linkcheck = true,
  linkcheck_ignore = [r"agupubs\.onlinelibrary\.wiley\.com", # This gives a 403 Forbidden
                      r"Decapodes\.jl/dev"], # 404, probably due to bad self-reference
  pages     = Any[
    "Decapodes.jl" => "index.md",
    "Overview" => "overview/overview.md",
    "Equations" => "equations/equations.md",
    "Vortices" => "navier_stokes/ns.md",
    "Harmonics" => "harmonics/harmonics.md",
    "Cahn-Hilliard" => "ch/cahn-hilliard.md",
    "Klausmeier" => "klausmeier/klausmeier.md",
    "CISM v2.1" => "cism/cism.md",
    "Glacial Flow" => "ice_dynamics/ice_dynamics.md",
    "Grigoriev Ice Cap" => "grigoriev/grigoriev.md", # Requires ice_dynamics
    "Budyko-Sellers-Halfar" => "bsh/budyko_sellers_halfar.md", # Requires ice_dynamics
    "Halfar-EBM-Water" => "ebm_melt/ebm_melt.md",
    "Halfar-NS" => "halmo/halmo.md", # Requires grigoriev
    "NHS" => "nhs/nhs_lite.md",
    "Pipe Flow" => "poiseuille/poiseuille.md",
    "Misc Features" => "bc/bc_debug.md", # Requires overview
    "ASCII Operators" => "ascii.md",
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
