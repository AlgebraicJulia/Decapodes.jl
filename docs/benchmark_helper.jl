"""
    insert_benchmarks!(file::String, benchmark_data::Dict)

Insert benchmark results into a markdown file at specified placeholder locations.
Placeholders should be in the format: <!-- BENCHMARK:benchmark_name -->

Parameters:
- file: Path to the markdown file
- benchmark_data: Dictionary mapping benchmark names to their result tables
"""
function insert_benchmarks!(file::String, benchmark_data::Dict{String,String})
	open(file, "w") do f
		foreach(eachline(f)) do line
			if startswith(line, "<!-- BENCHMARK:table -->")
				foreach(eachline(benchmark_data["table"])) do line
					println(f, line)
				end
			else
				println(f, line)
			end
		end 
	end
end
