srun --cpus-per-task=16 \
  --mem=16G \
  --time=2:00:00 \
  --output=.buildkite/build_%j.log \
  --unbuffered \
  .buildkite/jobscript.sh
