#!/bin/bash
# fine1 whole-layer launcher: the 98,280-voxel layer, arms CS, J5S, PS (and CF) one after another,
# each timed alone. --cap-h is required (1.5 x each arm's high projection, from the slice runs).
#   nohup caffeinate -i bash scripts/autoresearch/fine1/fine1_launch_whole.sh \
#       --cpp-bin <binary for CS> --cap-h CS=<h>,J5S=<h>,PS=<h> --resume > /dev/null 2>&1 &
#   bash scripts/autoresearch/fine1/fine1_launch_whole.sh --cpp-bin <binary> --cap-h CS=1,J5S=1,PS=1 --dry-run
# All options are those of fine1_launch.sh (read its header); the stage is fixed here. Nothing is
# parsed or scored until all end markers exist.
exec bash "$(dirname "$0")/fine1_launch.sh" --stage whole "$@"
