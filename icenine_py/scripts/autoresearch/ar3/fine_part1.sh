#!/bin/sh
# AR3 step 4, part 1: the arms side by side (5 heavy processes), each with a log and a .done marker.
# B015 and L run on the first 600 painted voxels of the voxel order only: with a box of 5/32 degree
# the fit's budget rule extends the budget many times (about 33,000 evaluations per fit in the
# pilot; median 22,145 in the run), so all painted voxels would cost hundreds of millions.
cd "$(dirname "$0")" || exit 1
unset VIRTUAL_ENV
L=../../../benchmarks/autoresearch/ar3/fine/cache/logs
mkdir -p "$L"
run() {
    name=$1
    shift
    (caffeinate -i uv run python -u run_fine.py "$@" > "$L/$name.log" 2>&1; touch "$L/$name.done") &
}
run B03 --arm B03
run B06 --arm B06
run C --arm C
run B015 --arm B015 --limit 600 --name B015
run L --arm L --limit 600 --name L
wait
