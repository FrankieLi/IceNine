#!/bin/bash
# tex1 whole-layer launcher: the 98,280-voxel layer, one arm per call (CSP, then RCP), each timed
# alone. --cap-h is required (plan: CSP=9, RCP=16 hours).
#   nohup caffeinate -i bash scripts/autoresearch/tex1/tex1_launch_whole.sh --arms CSP \
#       --cpp-bin <binary> --cap-h CSP=9 --allow-contended --continue-contended --resume \
#       > /dev/null 2>&1 &
#   nohup caffeinate -i bash scripts/autoresearch/tex1/tex1_launch_whole.sh --arms RCP \
#       --cap-h RCP=16 --allow-contended --continue-contended --resume > /dev/null 2>&1 &
#   bash scripts/autoresearch/tex1/tex1_launch_whole.sh --arms CSP --cpp-bin <binary> \
#       --cap-h CSP=9 --dry-run
# All options are those of tex1_launch.sh (read its header); the stage is fixed here. Nothing is
# parsed or scored until both end markers exist.
exec bash "$(dirname "$0")/tex1_launch.sh" --stage whole "$@"
