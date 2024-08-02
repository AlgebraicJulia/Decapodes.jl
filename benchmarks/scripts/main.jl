using DrWatson
@quickactivate :benchmarks

@info "Precompiling Julia"
using Pkg
Pkg.instantiate()
Pkg.precompile()
@info "Finished precompiling Julia"

using TOML

if length(ARGS) != 3
    error("Usage: 'sim_name' 'architecture' 'tag'")
end

# @info "Regenerating configuration files"
# run(`julia --threads=auto config_generate.jl`)

const sim_name = ARGS[1]
const arch = ARGS[2]
const tag = ARGS[3]
is_supported_arch(arch) || error("Architecture $(arch) is not in list $(supported_arches())")
sim_namedata = SimNameData(sim_name, arch, tag)
runname = tagfor_run(sim_namedata)

rm(resultsdir(sim_name), recursive=true, force=true)
mkpath(resultsdir(sim_name))

const config = config_path(sim_namedata)
dependency_ids = []

const cpu_job_args = "--partition=hpg-milan"
const cuda_job_args = "--partition=gpu --gres=gpu:a100:1"

if isfile(config)
    job_args = arch == "cpu" ? cpu_job_args : cuda_job_args
    config_data = TOML.parsefile(config)
    count = get_autoconfig_size(config_data)
    @info "Running $count $arch tasks"
    jobid = readchomp(`sbatch --array=1-$(count)%6 --parsable $(job_args) $(scriptsdir("array.sh")) $sim_name $arch $tag`)
    @info "Job ID for $(runname) is: $jobid"
    push!(dependency_ids, jobid)
else
    error("Configuration file for $(runname) could not be found")
end

run(`sbatch --dependency=afterok:$(join(dependency_ids, ",")) $(scriptsdir("final.sh")) $sim_name $arch $tag`)

