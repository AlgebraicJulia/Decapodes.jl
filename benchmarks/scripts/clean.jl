using DrWatson
@quickactivate "benchmarks"

# TODO: Can improve by seperating bench and final logs
file_regex = r"^.*log_(\d+)(?:_\d+)?\.txt$"

file_matches = filter(!isnothing, map(x -> match(file_regex, x), readdir(scriptsdir())))

for file in file_matches
  tgt_dir = scriptsdir("slurm_logs", "logs_"*file[1])
  mkpath(tgt_dir)
  mv(scriptsdir(file.match), joinpath(tgt_dir, file.match))
end
