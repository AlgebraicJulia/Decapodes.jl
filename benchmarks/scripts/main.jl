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

rm(resultsdir(sim_name), recursive=true, force=true)
mkpath(resultsdir(sim_name))

const cpu_config_file = get_config(sim_name, "cpu")
const cuda_config_file = get_config(sim_name, "cuda")

dependency_ids = []

# Run slurm scripts
if isfile(cpu_config_file)
    cpu_config = TOML.parsefile(cpu_config_file)
    cpu_count = length(cpu_config) - 1
    @info "Running $cpu_count CPU tasks"
    cpu_jobid = readchomp(`sbatch --array=1-$(cpu_count)%6 --parsable --partition=hpg-milan $(scriptsdir("array.sh")) $sim_name "cpu"`)
    @info "CPU Job ID is: $cpu_jobid"
    push!(dependency_ids, cpu_jobid)
end

if isfile(cuda_config_file)
    cuda_config = TOML.parsefile(cuda_config_file)
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

run(`julia --threads=auto clean.jl`)


