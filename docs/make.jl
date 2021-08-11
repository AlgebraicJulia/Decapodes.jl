using Documenter
using Literate

@info "Loading Decapods"
using Decapods
using Catlab
using Catlab.WiringDiagrams

# Set Literate.jl config if not being compiled on recognized service.
config = Dict{String,String}()
if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
  config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/Decapods.jl/blob/gh-pages/dev"
  config["repo_root_url"] = "https://github.com/AlgebraicJulia/Decapods.jl/blob/master/docs"
end

#const literate_dir = joinpath(@__DIR__, "..", "examples")
#const generated_dir = joinpath(@__DIR__, "src", "examples")
#
#for (root, dirs, files) in walkdir(literate_dir)
#  out_dir = joinpath(generated_dir, relpath(root, literate_dir))
#  for file in files
#    f,l = splitext(file)
#    if l == ".jl" && !startswith(f, "_")
#      Literate.markdown(joinpath(root, file), out_dir;
#        config=config, documenter=true, credit=false)
#      Literate.notebook(joinpath(root, file), out_dir;
#        execute=true, documenter=true, credit=false)
#    end
#  end
#end

@info "Building Documenter.jl docs"
makedocs(
  modules   = [Decapods],
  format    = Documenter.HTML(
    assets = ["assets/analytics.js"],
  ),
  sitename  = "Decapods.jl",
  doctest   = false,
  checkdocs = :none,
  pages     = Any[
    "Decapods.jl" => "index.md",
#    "Examples" => Any[
#      "examples/cfd_example.md"
#    ],
    "Library Reference" => "api.md"
  ]
)

@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/Decapods.jl.git",
  branch = "gh-pages"
)