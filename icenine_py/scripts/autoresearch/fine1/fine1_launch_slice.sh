#!/bin/bash
# fine1 slice launcher: the 2,000-voxel slice, arms CS, J5S, PS (and CF) one after another.
#   bash scripts/autoresearch/fine1/fine1_launch_slice.sh --cpp-bin <binary for CS> [--dry-run] \
#       [--resume] [--cf-bin <fixed binary>] [--cap-h CS=1,J5S=1,PS=1]
# All options are those of fine1_launch.sh (read its header); the stage is fixed here.
exec bash "$(dirname "$0")/fine1_launch.sh" --stage slice "$@"
