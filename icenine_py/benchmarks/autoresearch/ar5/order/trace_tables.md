<!-- table:order_noise -->
| sample | Y against A, after the pass (this step; order changed, draws aligned) | one re-draw of the kept pass tries against the full pass (step 7, variant 1) | Yl against Y, after the pass (this step; only the draws changed) |
|---|---|---|---|
| held | 1 / 6 | 4 / 10 |  |
| mid | 2 / 5 | 4 / 7 |  |
| grad | 0 / 0 | 15 / 12 |  |
| big | 3 / 0 | 2 / 6 | 5 / 12 |
<!-- /table:order_noise -->

<!-- table:order_trace_chain -->
| sample | seed steps at the same voxels from the start | first step whose taken end differs (counting from 1) | voxels that differ before the pass | voxels that differ after the pass | matched seed steps | of those, taken end differs | seed steps in one arm only (A / Y) | wrong ends, one-arm steps (A / Y) | wrong ends, all seed steps (A / Y) |
|---|---|---|---|---|---|---|---|---|---|
| held | 116 | 30 | 2,749 | 2,917 | 3,231 | 232 | 478 / 449 | 1 / 2 | 45 / 46 |
| mid | 65 | 52 | 3,049 | 3,329 | 1,173 | 115 | 434 / 437 | 4 / 2 | 15 / 13 |
| grad | 458 | 457 | 125 | 138 | 1,259 | 5 | 19 / 11 | 0 / 0 | 19 / 19 |
| big | 28 | 22 | 1,175 | 1,237 | 217 | 22 | 100 / 102 | 0 / 0 | 0 / 0 |
<!-- /table:order_trace_chain -->

<!-- table:order_trace_flips -->
| sample | voxel | Y against A | made | filled by, A / Y | error before the pass, A / Y (deg) | error after the pass, A / Y (deg) | seed steps at the voxel | region before the pass |
|---|---|---|---|---|---|---|---|---|
| held | 3084 | broken | in the pass | revisit / seed | 0.045 / 0.234 | 0.045 / 33.627 | Y-only seed step (no end reached 0.95): took hit 0.814 at 0.23 deg (highest hit 0.846 at 33.61 deg) | changed region with a matched seed step whose taken end differs (337 voxels) |
| held | 7689 | broken | before the pass | seed / revisit | 0.019 / 50.213 | 0.019 / 50.213 | A-only seed step (an end reached 0.95): took hit 0.996 at 0.02 deg | changed region with a matched seed step whose taken end differs (3 voxels) |
| held | 9669 | broken | before the pass | revisit / seed | 0.008 / 39.409 | 0.008 / 39.409 | Y-only seed step (no end reached 0.95): took hit 0.817 at 39.41 deg (highest hit 0.831 at 0.07 deg) | changed region with a matched seed step whose taken end differs (55 voxels) |
| held | 13777 | fixed | in the pass | revisit / seed | 0.284 / 0.004 | 31.315 / 0.004 | Y-only seed step (an end reached 0.95): took hit 0.989 at 0.00 deg | changed region with a matched seed step whose taken end differs (15 voxels) |
| held | 21444 | broken | before the pass | seed / revisit | 0.016 / 52.216 | 0.016 / 52.216 | A-only seed step (an end reached 0.95): took hit 0.990 at 0.02 deg | changed region with a matched seed step whose taken end differs (68 voxels) |
| held | 21987 | broken | in the pass | revisit / revisit | 43.482 / 43.482 | 0.03 / 43.482 | no seed step at the voxel | identical before the pass |
| held | 22475 | broken | before the pass | seed / revisit | 0.001 / 21.109 | 0.001 / 21.109 | A-only seed step (an end reached 0.95): took hit 0.991 at 0.00 deg | changed region with one-arm seed steps only (1 voxel) |
| mid | 1426 | broken | in the pass | neighbor / neighbor | 30.269 / 30.269 | 0.075 / 30.269 | no seed step at the voxel | identical before the pass |
| mid | 2078 | broken | before the pass | revisit / neighbor | 0.022 / 24.152 | 0.022 / 24.152 | no seed step at the voxel | changed region with a matched seed step whose taken end differs (836 voxels) |
| mid | 2541 | fixed | in the pass | revisit / revisit | 0.07 / 0.071 | 35.803 / 0.071 | no seed step at the voxel | changed region with one-arm seed steps only (24 voxels) |
| mid | 5209 | broken | in the pass | revisit / revisit | 0.06 / 0.098 | 0.06 / 53.321 | no seed step at the voxel | changed region with a matched seed step whose taken end differs (461 voxels) |
| mid | 8697 | fixed | in the pass | neighbor / revisit | 0.091 / 0.08 | 56.384 / 0.08 | no seed step at the voxel | changed region with a matched seed step whose taken end differs (256 voxels) |
| mid | 13013 | broken | before the pass | revisit / seed | 0.057 / 44.708 | 0.057 / 44.708 | Y-only seed step (no end reached 0.95): took hit 0.906 at 44.71 deg | changed region with one-arm seed steps only (6 voxels) |
| mid | 21630 | broken | in the pass | revisit / revisit | 0.056 / 0.046 | 0.056 / 54.741 | no seed step at the voxel | changed region with a matched seed step whose taken end differs (836 voxels) |
<!-- /table:order_trace_flips -->

<!-- table:order_trace_origin -->
| sample | flips after the pass, Y against A (fixed / broken) | made before the pass | made in the pass |
|---|---|---|---|
| held | 1 / 6 | 4 | 3 |
| mid | 2 / 5 | 2 | 5 |
| grad | 0 / 0 | 0 | 0 |
| big | 3 / 0 | 0 | 3 |
<!-- /table:order_trace_origin -->
