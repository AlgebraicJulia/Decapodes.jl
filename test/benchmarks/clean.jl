# TODO: Can improve by seperating bench and final logs
file_regex = r"^.*log_(\d+)(?:_\d+)?\.txt$"

file_matches = filter(!isnothing, map(x -> match(file_regex, x), readdir(".")))

for file in file_matches
  tgt_dir = joinpath("slurm_logs", "logs_"*file[1])
  mkpath(tgt_dir)
  mv(file.match, joinpath(tgt_dir, file.match))
end
