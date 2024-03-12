srun --cpus-per-task=16 \
  --mem=8G \
  --time=1:00:00 \
  --output=.buildkite/build_%j.log \
  --unbuffered \
  .buildkite/jobscript.sh
