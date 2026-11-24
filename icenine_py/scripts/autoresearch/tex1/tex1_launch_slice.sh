#!/bin/bash
# tex1 slice launcher: the 2,000-voxel slice, arm RCP (the pilot).
#   bash scripts/autoresearch/tex1/tex1_launch_slice.sh [--dry-run] [--resume] [--cap-h RCP=1] \
#       [--allow-contended --continue-contended]
# All options are those of tex1_launch.sh (read its header); the stage is fixed here.
exec bash "$(dirname "$0")/tex1_launch.sh" --stage slice "$@"
