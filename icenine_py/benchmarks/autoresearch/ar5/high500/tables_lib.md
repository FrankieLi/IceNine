<!-- table:high500_libcheck -->
| item (before the boundary pass) | library run LibJ5 | script arm J5l | equal |
|---|---|---|---|
| map: voxels that differ | 0 | - | yes |
| evaluations | 8,135,786 | 8,135,786 | yes |
| counters of the script run (36) | - | - | yes |
| per-voxel records (24,570) | - | - | yes |
| rescues: tried, rescued, tries | 11, 8, 27 | 11, 8, 27 | yes |
| MC refinements with the stall check (each: checks, stop, evaluations) | 1,689 | 1,689 | yes |
| stall checks; stopped; went on | 1,922; 402; 1,520 | 1,922; 402; 1,520 | yes |
| evaluations in stopped refinements | 466,575 | 466,575 | yes |

Equal before the pass: yes.
<!-- /table:high500_libcheck -->

<!-- table:high500_libpass -->
| item | library run LibJ5: pass in the run | script arm J5l: pass after the run |
|---|---|---|
| pass tries | 50,051 | 50,061 |
| pass evaluations | 3,053,111 | 3,053,721 |
| total with the pass | 11,188,897 | 11,189,507 |
| tries in the same (voxel, candidate) order from the start | 709 | - |
| map after the pass against the script's: voxels that differ in any bit / more than 0.001 degree apart / more than 1 degree apart | 2,157 / 2,132 / 10 | - |
| wrong voxels after the pass (more than 1 degree from the truth) | 39 | 37 |
| by voxel, script pass as the base: fixed / broken / wrong in both | 4 / 6 / 33 | - |
| median error after the pass, all 24,570 voxels (degrees) | 0.02913 | 0.02913 |
<!-- /table:high500_libpass -->
