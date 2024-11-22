"""
    insert_benchmarks!(file::String, benchmark_data::Dict)

Insert benchmark results into a markdown file at specified placeholder locations.
Placeholders should be in the format: <!-- BENCHMARK:benchmark_name -->

Parameters:
- file: Path to the markdown file
- benchmark_data: Dictionary mapping benchmark names to their result tables
"""
function insert_benchmarks!(file::String, benchmark_data::Dict{String,String})
	ks = keys(benchmark_data);
	open(file, "r+") do f
		foreach(readlines(file)) do line
			stub = filter(u -> startswith(line, "<!-- BENCHMARK:$u -->"), ks);
			if !isempty(stub)
				# get the most recent sims
				eligible_sims = readdir(benchmark_data[only(stub)], join=true)
				most_recent = filter(eligible_sims) do path
					last(splitpath(path)) == maximum(last.(splitpath.(eligible_sims)))
                end
				foreach(eachline(joinpath(only(most_recent), "default_output.md"))) do line
					println(f, line)
				end
			else
				println(f, line)
			end
		end 
	end
end
