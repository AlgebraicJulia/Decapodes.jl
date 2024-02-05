srun --cpus-per-task=32 \
  --mem=8G \
  --time=1:00:00 \
  --output=.buildkite/log_%jl.log \
  --unbuffered \
  .buildkite/jobscript.sh
