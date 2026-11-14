<!-- table:stallhit_libcheck -->
| item (before the boundary pass) | library run LibH | script arm Hl | equal |
|---|---|---|---|
| map: voxels that differ | 0 | - | yes |
| evaluations | 8,713,457 | 8,713,457 | yes |
| counters of the script run (36) | - | - | yes |
| per-voxel records (24,570) | - | - | yes |
| rescues: tried, rescued, tries | 12, 9, 29 | 12, 9, 29 | yes |
| MC refinements with the stall check (each: checks, stop, evaluations) | 1,597 | 1,597 | yes |
| stall checks; stopped; went on | 361; 197; 164 | 361; 197; 164 | yes |
| evaluations in stopped refinements | 748,128 | 748,128 | yes |

Equal before the pass: yes.
<!-- /table:stallhit_libcheck -->

<!-- table:stallhit_libpass -->
| item | library run LibH: pass in the run | script arm Hl: pass after the run |
|---|---|---|
| pass tries | 49,963 | 49,971 |
| pass evaluations | 3,047,743 | 3,048,231 |
| total with the pass | 11,761,200 | 11,761,688 |
| tries in the same (voxel, candidate) order from the start | 1,364 | - |
| map after the pass against the script's: voxels that differ in any bit / more than 0.001 degree apart / more than 1 degree apart | 2,109 / 2,068 / 10 | - |
| wrong voxels after the pass (more than 1 degree from the truth) | 36 | 38 |
| by voxel, script pass as the base: fixed / broken / wrong in both | 6 / 4 / 32 | - |
| median error after the pass, all 24,570 voxels (degrees) | 0.02913 | 0.02913 |
<!-- /table:stallhit_libpass -->
