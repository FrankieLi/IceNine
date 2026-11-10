<!-- table:passsets_l0 -->
| arm | edge_rules | pass_pairs | cpp_pass_pairs | fz_identical | pairs_pass_differs | overlap_equal | overlap_equal_pct |
|---|---|---|---|---|---|---|---|
| exact, one orientation per call (today) | off | 12,338 | 12,326 | 4,466 | 442 | 43,359 | 98.601 |
| exact, 32 FZ points per call (BatchedScreen) | off | 12,338 | 12,326 | 4,466 | 442 | 43,359 | 98.601 |
| interpolated, one FZ point per call (CppOmegaInterp) | off | 12,364 | 12,326 | 4,849 | 40 | 43,913 | 99.861 |
| interpolated, 32 FZ points per call (both) | off | 12,364 | 12,326 | 4,849 | 40 | 43,912 | 99.859 |
| exact, one orientation per call (today) | on | 12,298 | 12,326 | 4,500 | 402 | 43,416 | 98.731 |
| exact, 32 FZ points per call (BatchedScreen) | on | 12,298 | 12,326 | 4,500 | 402 | 43,416 | 98.731 |
| interpolated, one FZ point per call (CppOmegaInterp) | on | 12,325 | 12,326 | 4,885 | 1 | 43,970 | 99.991 |
| interpolated, 32 FZ points per call (both) | on | 12,325 | 12,326 | 4,885 | 1 | 43,969 | 99.989 |
<!-- /table:passsets_l0 -->

<!-- table:speed_screen -->
| arm | us_total | us_solve | us_overlap_c | us_rest | speedup |
|---|---|---|---|---|---|
| exact, one orientation per call (today) | 388.4 | 61.7 | 8.6 | 333.2 | 1.0 |
| exact, 32 FZ points per call (BatchedScreen) | 22.2 | 0.6 | 7.7 | 14.1 | 17.5 |
| interpolated, one FZ point per call (CppOmegaInterp) | 80.5 | 6.8 | 7.4 | 67.2 | 4.8 |
| interpolated, 32 FZ points per call (both) | 22.7 | 0.3 | 7.7 | 15.0 | 17.1 |
<!-- /table:speed_screen -->

<!-- table:speed_voxel -->
| switches | n | l0_screen_us | l0_screen_s | levels012_s | wall_s | evals | contended |
|---|---|---|---|---|---|---|---|
| off | 6 | 415.2 | 18.3 | 115.6 | 119.3 | 225,883 | 5 |
| BatchedScreen | 6 | 24.0 | 1.1 | 89.6 | 92.6 | 225,883 | 4 |
| CppOmegaInterp | 6 | 88.3 | 3.9 | 94.2 | 97.5 | 225,648 | 3 |
| both | 6 | 24.7 | 1.1 | 88.9 | 92.0 | 225,648 | 1 |
<!-- /table:speed_voxel -->

<!-- table:speed_voxel_quiet -->
| arm | voxel | label | wall_s | evals | l0_screen_us | us_per_eval | load |
|---|---|---|---|---|---|---|---|
| cpp_edge | 15901 | single-worker | 122.5 | 237,890 | 406.9 | 493 | 1.39 |
| cpp_edge | 5242 | single-worker | 100.8 | 199,386 | 398.7 | 486 | 1.29 |
| cpp_edge_interp_batched | 15901 | single-worker | 96.2 | 238,868 | 24.2 | 381 | 1.35 |
| cpp_edge_interp_batched | 5242 | single-worker | 75.6 | 199,195 | 24.1 | 360 | 1.19 |
<!-- /table:speed_voxel_quiet -->
