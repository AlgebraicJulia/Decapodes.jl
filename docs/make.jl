using Documenter
using Literate

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

const literate_dir = joinpath(@__DIR__, "..", "examples")
const generated_dir = joinpath(@__DIR__, "src", "examples")

for (root, dirs, files) in walkdir(literate_dir)
  out_dir = joinpath(generated_dir, relpath(root, literate_dir))
  for file in files
    f,l = splitext(file)
    if l == ".jl" && !startswith(f, "_")
      Literate.markdown(joinpath(root, file), out_dir;
        config=config, documenter=true, credit=false)
      Literate.notebook(joinpath(root, file), out_dir;
        execute=true, documenter=true, credit=false)
    end
  end
end

pages = Pair{String, String}[]
dirs = Dict("physics"    => "Physics"
            ,"chemistry" => "Chemistry"
            ,"biology"   => "Biology")
push!(f, "Library Reference" => "api.md")
for d in keys(dirs), f in cd(readdir, d)
    push!(pages, f => dirs[d] * "/" * readlines(f)[1])
end

@info "Building Documenter.jl docs"
makedocs(
  modules   = [Decapodes],
  format    = Documenter.HTML(
    assets = ["assets/analytics.js"],
  ),
  sitename  = "Decapodes.jl",
  doctest   = false,
  checkdocs = :none,
  pages = pages
 )

@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/Decapodes.jl.git",
  branch = "gh-pages",
  devbranch = "main"
)
