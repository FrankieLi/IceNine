#!/bin/bash
# Centroid task step 3 (clean images only). Run from icenine_py/. All wall times are contended (10 workers).
cd "$(dirname "$0")/../../.."
L=benchmarks/phase_d_seed_diag/centroid/cache/logs
export RAY_VARIANTS=clean
date +%T > $L/start
# (a) step 1 peak level, 30 voxels, clean
RAY_PEAK_DIR=peak_centroid uv run python scripts/phase_d/ray_overlap_study.py run --step 1 --workers 10 > $L/a_peak.log 2>&1
touch $L/a_peak.done
# (b) 30-voxel coarse test, 4 arms, clean
RAY_COARSE_DIR=coarse30_centroid uv run python scripts/phase_d/ray_overlap_study.py run --step 2 --workers 10 --arms old,cpp,f0_R2,f0_A > $L/b_coarse.log 2>&1
touch $L/b_coarse.done
# (c) E0 200 voxels clean, min_sin_eta 0, arms old (off), cpp (on), R2, A
for arm in off on rayR2 rayA; do
  uv run python scripts/phase_d/e0_cpprule.py run --arm $arm --msin 0 --workers 10 --tag centroid --variant clean > $L/c_e0_$arm.log 2>&1
done
touch $L/c_e0.done
# (d) six full-sample voxels, old arm (6 single-process runs in parallel, contended)
for v in 15901 17034 18558 19889 22867 5242; do
  uv run python scripts/phase_d/seed_diag.py run --voxel $v --images full --tag base_centroid > $L/d_six_$v.log 2>&1 &
done
wait
touch $L/d_six.done
date +%T > $L/end
touch $L/all.done
