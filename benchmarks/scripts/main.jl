using DrWatson
@quickactivate :benchmarks

@info "Precompiling Julia"
using Pkg
Pkg.instantiate()
Pkg.precompile()
@info "Finished precompiling Julia"

using TOML

if length(ARGS) != 1
    error("Usage: 'sim_name'")
end

const sim_name = ARGS[1]
simdir(args...) = srcdir(sim_name, args...)

# const tracker = simdir("clean.txt")

mkpath(paramsdir(sim_name))
rm(resultsdir(sim_name), recursive=true, force=true)
mkpath(resultsdir(sim_name))

# Force parameter regeneration
# if !isfile(tracker)
#     for file in readdir(paramsdir(), join=true)
#         rm(file)
#     end
#     open(tracker, "w") do file
#         write(file, "Configuration has not been updated, tuning parameters still valid")
#     end
# end

const cpu_config_name = get_configname(sim_name, "cpu")
const cuda_config_name = get_configname(sim_name, "cuda")

dependency_ids = []

# Run slurm scripts
if isfile(simdir(cpu_config_name))
    cpu_config = TOML.parsefile(simdir(cpu_config_name))
    cpu_count = length(cpu_config) - 1
    @info "Running $cpu_count CPU tasks"
    cpu_jobid = readchomp(`sbatch --array=1-$(cpu_count)%6 --parsable --partition=hpg-milan $(scriptsdir("array.sh")) $sim_name "cpu"`)
    @info "CPU Job ID is: $cpu_jobid"
    push!(dependency_ids, cpu_jobid)
end

if isfile(simdir(cuda_config_name))
    cuda_config = TOML.parsefile(simdir(cuda_config_name))
    cuda_count = length(cuda_config) - 1
    @info "Running $cuda_count CUDA tasks"
    cuda_jobid=readchomp(`sbatch --array=1-$(cuda_count)%6 --parsable --partition=gpu --gres=gpu:a100:1 $(scriptsdir("array.sh")) $sim_name "cuda"`)
    @info "CUDA Job ID is: $cuda_jobid"
    push!(dependency_ids, cuda_jobid)
end

if isempty(dependency_ids)
    error("No configuration information found for $sim_name")
end

run(`sbatch --dependency=afterok:$(join(dependency_ids, ",")) $(scriptsdir("final.sh")) $sim_name`)




