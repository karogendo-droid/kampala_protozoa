#!/bin/bash
# One-command launcher. Does not remove, move, or overwrite any sequencing data.
set -euo pipefail
PROJECT="$HOME/bioinformatics/kampala_protozoa"
if [ ! -f "$PROJECT/scripts/10_batch_pipeline.py" ]; then
  echo "STOP: Could not find the batch pipeline in $PROJECT/scripts"
  exit 2
fi
cd "$PROJECT"
python3 scripts/10_batch_pipeline.py "$@"
