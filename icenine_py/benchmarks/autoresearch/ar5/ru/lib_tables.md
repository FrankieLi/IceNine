<!-- table:ru_lib -->
| L1: compared before the boundary pass | script arm RUl | library, key on | verdict |
|---|---|---|---|
| map (decompressed .mic bytes) | - | - | equal |
| evaluations | 4,721,521 | 4,721,521 | equal |
| counters of the script run | 36 | all | equal |
| per-voxel records | 24570 | 24570 | equal |
| seed steps (try rows) | 1692 | 1692 | equal |
| seed steps with a try | 1691 | 1691 | equal |
| tries that succeeded | 1379 | 1379 | equal |
| steps with the PXF step | 313 | 313 | equal |
| MC refinements with checks | 470 | 470 | equal |
<!-- /table:ru_lib -->

<!-- table:ru_lib_l2 -->
| L2: library form, mid | LibPXF | LibRU | verdict |
|---|---|---|---|
| evaluations before the pass | 6,731,510 | 4,721,521 | 29.9% fewer (met) |
| wrong voxels after the pass | 31 | 32 | broken 11, fixed 10, net 1 (met) |
| median error after the pass (deg) | 0.028710 | 0.028223 | -0.000488 (met) |
<!-- /table:ru_lib_l2 -->
