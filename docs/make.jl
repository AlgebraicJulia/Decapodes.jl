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

const literate_dir = joinpath(@__DIR__, "..", "examples")
const generated_dir = joinpath(@__DIR__, "src", "examples")

@info "Building literate files"
for (root, dirs, files) in walkdir(literate_dir)
  out_dir = joinpath(generated_dir, relpath(root, literate_dir))
  pmap(files) do file
    f,l = splitext(file)
    if l == ".jl" && !startswith(f, "_")
      Literate.markdown(joinpath(root, file), out_dir;
        config=config, documenter=true, credit=false)
      Literate.notebook(joinpath(root, file), out_dir;
        execute=true, documenter=true, credit=false)
    end
  end
end
@info "Completed literate"

pages = Any[]
push!(pages, "Decapodes.jl"      => "index.md")
push!(pages, "Overview"          => "overview.md")
push!(pages, "Equations"         => "equations.md")
push!(pages, "BC Debug"          => "bc_debug.md")
push!(pages, "ASCII Operators"   => "ascii.md")
dirs = Dict(
	   "physics"  => "Physics"
	   ,"biology"  => "Biology"
	   "climate"  => "Climate")
for d in keys(dirs)
  dir   = joinpath(@__DIR__, "src", d)
  files = readdir(dir)
  push!(pages, dirs[d] => joinpath.(d, files))
end

push!(pages, "Examples" => [
	"Brusselator"        => "examples/chemistry/brusselator.md",
	"Brusselator Teapot" => "examples/chemistry/brusselator_teapot.md",
	"Gray-Scott"         => "examples/chemistry/gray_scott.md",
	"Budyko-Sellers"     => "examples/climate/budyko_sellers.md",
	"Budyko-Sellers-Halfar" => "examples/climate/budyko_sellers_halfar.md",
	"Burger" => "examples/diff_adv/burger.md",
])
push!(pages, "Canonical Models"  => "canon.md")
push!(pages, "Library Reference" => "api.md")

@info "Building Documenter.jl docs"
makedocs(
  modules   = [Decapodes],
  format    = Documenter.HTML(
    assets = ["assets/analytics.js"],
  ),
  remotes   = nothing,
  sitename  = "Decapodes.jl",
  doctest   = false,
  checkdocs = :none,
  pages     = pages)


@info "Deploying docs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/Decapodes.jl.git",
  branch = "gh-pages",
  devbranch = "main"
)
