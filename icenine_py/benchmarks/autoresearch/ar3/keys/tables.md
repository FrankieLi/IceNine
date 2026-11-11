<!-- table:keys_evals -->
| run | code | seeds | reseeds | grid | seed_step | grid_fallback | neighbour | revisit | pass | total | wrong_before | wrong_after |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D4 | script | 38 | 291 | 8 | 4016458 | 2452689 | 3479526 | 863042 | 1131977 | 11943692 | 68 | 21 |
| D4 | library | 38 | 291 | 8 | 4016458 | 2452689 | 3479526 | 863042 | 1131977 | 11943692 | 68 | 22 |
| S5 | script | 35 | 314 | 0 | 7469574 | 0 | 3513895 | 884040 | 1132343 | 12999852 | 71 | 30 |
| S5 | library | 35 | 314 | 0 | 7469574 | 0 | 3513895 | 884040 | 1132343 | 12999852 | 71 | 30 |
| S8 | script | 37 | 292 | 0 | 2856712 | 0 | 3491437 | 853714 | 1138565 | 8340428 | 67 | 24 |
| S8 | library | 37 | 292 | 0 | 2856712 | 0 | 3491437 | 853714 | 1138504 | 8340367 | 67 | 24 |
<!-- /table:keys_evals -->

<!-- table:keys_exact -->
| run | number | library | stored | verdict |
|---|---|---|---|---|
| D4L | seeds | 38 | 38 | exact |
| D4L | reseeds | 291 | 291 | exact |
| D4L | grid | 8 | 8 | exact |
| D4L | evals_before_pass | 10,811,715 | 10,811,715 | exact |
| D4L | counters (all numeric, no wall times) | 0 differ | 36 compared | exact |
| D4L | per-voxel records | equal | 24,570 voxels | exact |
| D4L | map before the pass vs stored recon.mic | equal | 6,535,568 bytes | exact |
| D4L | wrapper: fits | 32,562 | 32,562 | exact |
| D4L | wrapper: mc_fallbacks | 3,194 | 3,194 | exact |
| D4L | wrapper: fallback_won | 4 | 4 | exact |
| S5L | seeds | 35 | 35 | exact |
| S5L | reseeds | 314 | 314 | exact |
| S5L | grid | 0 | 0 | exact |
| S5L | evals_before_pass | 11,867,509 | 11,867,509 | exact |
| S5L | counters (all numeric, no wall times) | 0 differ | 33 compared | exact |
| S5L | per-voxel records | equal | 24,570 voxels | exact |
| S5L | map before the pass vs stored recon.mic | equal | 6,535,568 bytes | exact |
| S5L | wrapper: fits | 32,673 | 32,673 | exact |
| S5L | wrapper: mc_fallbacks | 3,235 | 3,235 | exact |
| S5L | wrapper: fallback_won | 3 | 3 | exact |
| S8L | seeds | 37 | 37 | exact |
| S8L | reseeds | 292 | 292 | exact |
| S8L | grid | 0 | 0 | exact |
| S8L | evals_before_pass | 7,201,863 | 7,201,863 | exact |
| S8L | counters (all numeric, no wall times) | 0 differ | 33 compared | exact |
| S8L | per-voxel records | equal | 24,570 voxels | exact |
| S8L | map before the pass vs stored recon.mic | equal | 6,535,568 bytes | exact |
| S8L | early stop: stopped_early | 322 | 322 | exact |
| S8L | early stop: seed_steps | 329 | 329 | exact |
| S8L | early stop: seed_step_evals | 2,856,712 | 2,856,712 | exact |
<!-- /table:keys_exact -->

<!-- table:keys_fresh -->
| run | compared | tries_a | tries_b | evals_a | evals_b | voxels_differing | verdict |
|---|---|---|---|---|---|---|---|
| D4L | (i) rerun's in-run pass vs step 10's in-run pass | 18,557 | 18,557 | 1,131,977 | 1,131,977 | 0 | equal |
| D4L | (iii) same process, rounded map vs script tries | 18,557 | 18,557 | 1,131,977 | 1,131,977 | 0 | equal |
| D4L | (iv) fresh process, unrounded map vs in-run pass | 18,557 | 18,557 | 1,131,977 | 1,131,977 | 0 | equal |
| D4L | (v) same process, rounded map vs in-run pass | 18,557 | 18,557 | 1,131,977 | 1,131,977 | 11087 | differ |
| D4L | (ii) before-pass .mic vs step 10's (bytes) |  |  |  |  |  | equal |
| D4L | (ii) map read back vs the script's before_*.npy |  |  |  |  |  | equal |
| S8L | (i) rerun's in-run pass vs step 10's in-run pass | 18,664 | 18,664 | 1,138,504 | 1,138,504 | 0 | equal |
| S8L | (iii) same process, rounded map vs script tries | 18,665 | 18,665 | 1,138,565 | 1,138,565 | 0 | equal |
| S8L | (iv) fresh process, unrounded map vs in-run pass | 18,664 | 18,664 | 1,138,504 | 1,138,504 | 0 | equal |
| S8L | (v) same process, rounded map vs in-run pass | 18,665 | 18,664 | 1,138,565 | 1,138,504 | 11054 | differ |
| S8L | (ii) before-pass .mic vs step 10's (bytes) |  |  |  |  |  | equal |
| S8L | (ii) map read back vs the script's before_*.npy |  |  |  |  |  | equal |
<!-- /table:keys_fresh -->

<!-- table:keys_measured -->
| sample | run | evals_with_pass | wrong_before | before_ci_pct | wrong_after | after_ci_pct |
|---|---|---|---|---|---|---|
| big | D4 | 11943692 | 68 | 0.22 to 0.35 | 21 | 0.06 to 0.13 |
| big | S5 | 12999852 | 71 | 0.23 to 0.36 | 30 | 0.09 to 0.17 |
| big | S8 | 8340428 | 67 | 0.21 to 0.35 | 24 | 0.07 to 0.15 |
| grad | D4 | 38807186 | 225 | 0.80 to 1.04 | 74 | 0.24 to 0.38 |
| grad | S5 | 20629760 | 203 | 0.72 to 0.95 | 62 | 0.20 to 0.32 |
| grad | S8 | 18000230 | 210 | 0.75 to 0.98 | 63 | 0.20 to 0.33 |
| mid | D4 | 40560874 | 94 | 0.31 to 0.47 | 48 | 0.15 to 0.26 |
| mid | S5 | 39093613 | 87 | 0.29 to 0.44 | 34 | 0.10 to 0.19 |
| mid | S8 | 26443074 | 84 | 0.28 to 0.42 | 31 | 0.09 to 0.18 |
<!-- /table:keys_measured -->

<!-- table:keys_off -->
| file | parent vs branch | bytes |
|---|---|---|
| recon.mic | equal | 530,981 |
| orient.npy | equal | 144,128 |
| records.json | equal | 505,215 |
| stats.json | equal | 843 |
| logs.json | equal | 31,917 |
<!-- /table:keys_off -->

<!-- table:keys_pass1 -->
| stored map | tries | tries_stored | evals | evals_stored | tries_differing | voxels_changed | other_end |
|---|---|---|---|---|---|---|---|
| D4 | 18,557 | 18,557 | 1,131,977 | 1,131,977 | 0 | 1057 | 0 |
| S5 | 18,563 | 18,563 | 1,132,343 | 1,132,343 | 0 | 1038 | 0 |
| S8 | 18,665 | 18,665 | 1,138,565 | 1,138,565 | 0 | 1056 | 0 |
<!-- /table:keys_pass1 -->

<!-- table:keys_pass2 -->
| run | tries | tries_stored | evals | evals_stored | same_voxel_and_cand | same_end | wrong_after | wrong_after_stored | fixed | broken | both_wrong |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D4L | 18,557 | 18,557 | 1,131,977 | 1,131,977 | 18556 | 14 | 22 | 21 | 0 | 1 | 21 |
| S5L | 18,563 | 18,563 | 1,132,343 | 1,132,343 | 18563 | 11 | 30 | 30 | 0 | 0 | 30 |
| S8L | 18,664 | 18,665 | 1,138,504 | 1,138,565 | 18661 | 5 | 24 | 24 | 5 | 5 | 19 |
<!-- /table:keys_pass2 -->

<!-- table:keys_uncounted -->
| run | boundary_score_evals | early_stop_hit_evals | wall_s | pass_wall_s |
|---|---|---|---|---|
| D4L | 11,087 | 0 | 1,090 | 125.4 |
| S5L | 11,049 | 0 | 1,156 | 123.5 |
| S8L | 11,054 | 465 | 865 | 128.0 |
<!-- /table:keys_uncounted -->
