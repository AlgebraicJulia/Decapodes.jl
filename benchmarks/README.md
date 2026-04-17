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

   ```julia
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

which auto-activate the project, enable local path handling from DrWatson and provide several helper functions. Note that `@quickactivate :benchmarks` is only usable from within the `benchmarks` directory.

## Establishing Simulation Configurations

To establish a set of configurations for a simulation, you should create a file `src/$physics/config.toml` with the below structure.

1. For a given physics (heat, brusselator, etc), entries are structured as `[$physics.$architecture.$tag]`, e.g. `[heat.cpu.example]`.
2. Under an entry, list all the parameters desired. This should be structured as either `$param_name = [val1, val2, ...]` or `$param_name = val`.
3. You should always include a `code_target`, which takes either `CPUTarget` or `CUDATarget` as a string.

## Creating a Simulation File

These benchmarks depend upon you to create the simulation files to be benchmarked. For a certain simulation named ```example```, the simulation file would be ```src/example/example.jl```.

Always start the file with the following, as it provides you access to helper functions located in ```src/helpers```.

```julia
using DrWatson
@quickactivate :benchmarks
```

Additionally, this file should contain the following functions.

1. ```setup_config```, which will take in a ```Dict{String, Any}``` with the provided parameter names pointing to that task's provided configuration values. It's up to you to take these values, process them and then organize them to be passed along.
2. ```setup_benchmark```, which will create the Decapode and run ```eval(gensim())``` on it and return that evaluated function.
3. ```create_mesh```, which will create the mesh upon which the simulation will run and also initialize the initial conditions and any constants/parameters. Return the mesh, initial conditions and constants/parameters in that order.
4. ```create_simulate```, which will take the generated mesh and evaluated function and run the simulate function. Return the resulting function.
5. ```run_simulation```, which will take the resulting simulation function, initial conditions and constants/parameters and run the solve. Return the result of the solve.

## Running the Benchmarks

Use the `main_config.toml` to list out which physics configuration entries you would like to be run. These entries correspond one-to-one with the entries in the physics configurations. However, these `main_config.toml` entries can take optional arguements to customize how the simulations are run. Currently supported arguments are:

1. `slurm_args`, passed as a vector of strings for each seperate `sbatch` argument. These arguments will be added to each job for that entry.
2. `concur_jobs`, passed as an integer that determines how many jobs can run at a time for that configuration.

Once done, simply run `scripts/main.sh`. 

As another option, the name of a specific configuration entry in the `main_config.toml` can be passed to `main.sh` to only run that job, e.g. `main.sh heat cpu example`.

**Warning**: The caller of this script should be in the `benchmarks` directory.

## Data Collection and Processing

Once a simulation is run, their outputs will be saved in `data/sims/$physics`. The benchmark JSON files will contain the full result of the benchmarking run while the stat JLD2 files will contain the solve result statistics from DEStats from OrdinaryDiffEq.jl.

These files will then be processed and have their data stored in `data/exp_pro/$physics/$slurm_job_id/autogen`. Each result file in that directory will contain the processed results from one single task.

These result files can then be collected with `collect_results` on the `autogen` directory and post-processed with `DataFrames.jl` functions to create the desired tables. An example of this kind of post-processing script is included in `scripts/post_processing`.

Because these files are expected to be collected using DrWatson's `collect_results` , there is no intended correlation between a task and its result file name.

**Warning**: Not all information from the benchmarking run is saved to the result files and any files in `data/sims/$physics`, specifically the JSON and JLD2 files mentioned above, will be deleted upon the next benchmark run for that physics. On the other hand, result files in the `autogen` directory mentioned before will never be deleted by the benchmarking.

## Testing

The benchmarks have a few tests that can be used to establish the robustness of the system. To run them, activate the `benchmarks` environment and then enter `test`.
