# Benchmarks

## DrWatson Initialization

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> benchmarks

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate :benchmarks
```
which auto-activate the project, enable local path handling from DrWatson and provide several helper functions.

## Establishing Configurations

To establish a set of configurations for a simulation, you should use ```config_generate.jl```, located in the ```scripts``` directory. This will read either the included ```main_config.toml``` or a user provided file, when the file name is provided as a command line argument, and distribute the configuration information to multiple tasks. 

To use a custom TOML, please follow all the steps below:

0. The user provided TOML must located in the ```src``` directory.
1. For a given simulation (heat, brusselator, etc.) and intended architeture (cpu, cuda), include a group in the TOML as `[sim_name.arch]`, where `sim_name` and `arch` are respective names from the previous groups.
2. Under a group, list all the parameters desired. This should be structured as either `param_name = [val1, val2, val3...]` or `param_name = val`.
3. Highly recommended parameters to include are ```code_target```, which takes either ```CPUTarget``` or ```CUDATarget```.

```config_generate.jl``` will generate the required directories automatically. If you were to include ```[heat.cpu]``` in your TOML, the resulting configuration would be in ```src/heat/heat_cpu.toml```.

Please view ```main_config.toml``` as a guiding example on how to craft your own TOML. 

**Warning**: ```config_generate.jl``` is not called automatically so it is up to you to run the script before launching benchmarks. 

**Warning**: Do not call ```config_generate.jl``` for an actively running simulation since this may unexpectedly affect the running process.


## Creating a Simulation File

These benchmarks depend upon you to create the simulation files to be benchmarked. For a certain simulation named ```example```, include a file in ```src/example/example.jl```. Always start the file with the following, as it provides you access to helper functions located in ```src/helpers```.
```julia
using DrWatson
@quickactivate :benchmarks
```

Additionally, this file should contain the following functions.

0. ```setup_config```, which will take in a ```Dict{String, Any}``` with the provided parameter names pointing to that task's provided configuration values. It is up to you to take these values, process them and then organization them to be passed along.
1. ```setup_benchmark```, which will create the Decapode and run ```eval(gensim())``` on it. Return the evaluated function.
2. ```create_mesh```, which will create the mesh upon which the simulation will run and also initialize the initial conditions and any constants/parameters. Return the mesh, initial conditions and constants/parameters in that order.
3. ```create_simulate```, which will take the generated mesh and evaluated function and run the simulate function. Return the resulting function.
4. ```run_simulation```, which will take the resulting simulation function, initial conditions and constants/parameters and run the solve. Return the result of the solve.

## Running the Benchmarks

Simply run ```scripts/main.sh``` and provide it the name of the desired simulation to benchmark. Note this must match the name provided in the user's config.

**Warning**: The caller of this script should be in the `benchmarks` directory.

Note that it may be the case that a benchmark run will hit an Out-of-Memory error. If that occurs, then you may update the reserved memory in ```scripts/array.sh```. Similarly, you may increase the time requested to run the benchmark.

## Data Collection and Processing

Once a simulation is run, their outputs will be saved in ```data/sims/"sim_name"```. The benchmark JSON files will contain the full result of the benchmarking run while the stat JLD2 files will contain the solve result statistics from DEStats.

These files will then be processed and have aggregated data stored in ```data/exp_pro/"sim_name"/"slurm_job_id"/autogen```. Each result file in that directory will have the entire results from one single task. Because these files are expected to be collected using DrWatson's ```collect_results``` on the ```autogen``` directory, there is no intended correlation between a task and its result file. 

**Warning**: Note that not all information from the benchmarking run is saved to the result files and any files in ```data/sims/"sim_name"``` will be deleted upon the next benchmark run. On the other hand, result files in the ```autogen``` directory mentioned before will never be deleted by the benchmarking.

An example Markdown file is output in ```data/exp_pro/"sim_name"/"slurm_job_id"``` for user inspection, called ```final.md```.

## Testing

The benchmarks are have a few tests that will eventually be used to establish the robustness of the system. These can be run by setting the active project to ```benchmarks```, either through the package manager directory or DrWatson, and then running ```test```.