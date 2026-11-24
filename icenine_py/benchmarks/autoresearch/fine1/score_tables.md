Reading line: one row per arm on the whole fine layer (98,280 voxels, clean images). Wrong is over 1 degree from the truth, cubic symmetry; a voxel with no output counts as wrong. All times contended.

| arm | workers | evaluations | wrong | share wrong (95% Wilson) | missing | median error (deg) | wall (s) |
|---|---|---|---|---|---|---|---|
| CS | 1 | 1,670,595,806 | 15604 | 15.877% (15.650, 16.107) | 10845 | 0.0852 | 91206 |
| CF | 1 | 1,362,495,605 | 12777 | 13.001% (12.792, 13.212) | 9926 | 0.0753 | 71996 |
| CSP | 12 | 1,380,035,950 | 14914 | 15.175% (14.952, 15.401) | 10411 | 0.0687 | 10267 |
| CFP | 12 | 1,370,501,354 | 13920 | 14.164% (13.947, 14.383) | 10590 | 0.0708 | 10213 |
| J5S | 1 | 32,579,814 | 561 | 0.571% (0.526, 0.620) | 0 | 0.0295 | 7557 |
| PS | 1 | 2,832,907,439 | 714 | 0.726% (0.675, 0.782) | 0 | 0.0617 | 183303 |

Reading line: the same count on the 85676 voxels that every arm has an output for.

| arm | wrong on the common voxels | share (95% Wilson) |
|---|---|---|
| CS | 4681 | 5.464% (5.313, 5.618) |
| CF | 2734 | 3.191% (3.075, 3.311) |
| CSP | 4416 | 5.154% (5.008, 5.304) |
| CFP | 3241 | 3.783% (3.657, 3.913) |
| J5S | 453 | 0.529% (0.482, 0.580) |
| PS | 616 | 0.719% (0.665, 0.778) |

Reading line: voxel-by-voxel comparison of two maps; `differ` counts voxels whose two orientations are more than 1e-6 degree apart.

| pair | identical | differ | differ over 0.05 deg | differ over 1 deg | wrong in both | only in first | only in second |
|---|---|---|---|---|---|---|---|
| CS vs CSP | 29 | 98251 | 60302 | 15374 | 12561 | 3043 | 2353 |
| CF vs CFP | 27 | 98253 | 56917 | 14946 | 10215 | 2562 | 3705 |
| CS vs CF | 2065 | 96215 | 55435 | 15771 | 10834 | 4770 | 1943 |
| CS vs PS | 0 | 98280 | 76886 | 15454 | 483 | 15121 | 231 |
| CS vs J5S | 0 | 98280 | 72238 | 15455 | 406 | 15198 | 155 |
| PS vs J5S | 0 | 98280 | 61855 | 744 | 266 | 448 | 295 |
| CSP vs CFP | 4134 | 94146 | 54934 | 14995 | 11879 | 3035 | 2041 |
