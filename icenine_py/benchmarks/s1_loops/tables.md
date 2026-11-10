<!-- table:s1_batched -->
| check | modes | (voxel, mode, variant) rows | identical |
|---|---|---|---|
| BatchedCostEvaluator (CPU compact and fixed) against the base evaluate | cpp_rule, edge, edge_cpp_rule, local, min_sin_eta | 40 | 40 |
<!-- /table:s1_batched -->

<!-- table:s1_coverage -->
| voxel | cost | omega rows | in [low - width, low) | above the range | further below | in gap bins | frames used |
|---|---|---|---|---|---|---|---|
| three0 | local | 202504 | 556 | 50302 | 50394 | 0 | 180 |
| three0 | screen | 45556 | 67 | 11283 | 11428 | 0 | 180 |
| full15901 | local | 206636 | 391 | 52919 | 50008 | 0 | 180 |
| full15901 | screen | 48492 | 64 | 13075 | 11107 | 0 | 180 |
| full17034 | local | 209096 | 769 | 52597 | 51182 | 0 | 180 |
| full17034 | screen | 49888 | 287 | 12570 | 12087 | 0 | 180 |
| full18558 | local | 207656 | 470 | 51371 | 51987 | 0 | 180 |
| full18558 | screen | 48332 | 122 | 12429 | 11615 | 0 | 180 |
<!-- /table:s1_coverage -->

<!-- table:s1_e0 -->
| run | voxels | wrong (> 1 deg) | wrong, the 200 voxels of e0_s3off | evaluations, global | evaluations, local |
|---|---|---|---|---|---|
| base commit | 208 | 37 | 35 | 9256059 | 858079 |
| S1 | 208 | 37 | 35 | 9256059 | 858079 |
<!-- /table:s1_e0 -->

<!-- table:s1_parity -->
| mode | orientations (4 voxels) | field arrays | differing values |
|---|---|---|---|
| cpp_rule | 4000 | 32 | 0 |
| edge | 4000 | 32 | 0 |
| edge_cpp_rule | 4000 | 32 | 0 |
| fallback_cpp_edge | 200 | 32 | 0 |
| fallback_local | 200 | 32 | 0 |
| fallback_screen | 200 | 32 | 0 |
| local | 4000 | 32 | 0 |
| min_sin_eta | 4000 | 32 | 0 |
| ray_c_local | 4000 | 32 | 0 |
| ray_c_screen | 4000 | 32 | 0 |
| ray_py_screen | 4000 | 32 | 0 |
| screen_center | 4000 | 32 | 0 |
| screen_edge | 4000 | 32 | 0 |
| screen_vertex0 | 4000 | 32 | 0 |
| soft | 200 | 32 | 0 |
<!-- /table:s1_parity -->

<!-- table:s1_per_call -->
| arm | base commit (us) | S1 (us) | saved (us) | base / S1 | base min | S1 min | rounds | aten ops per call, base / S1 |
|---|---|---|---|---|---|---|---|---|
| screen (pixel radius 3, Q-max 5, C++ rule) | 244.8 | 215.0 | 29.8 | 1.14 | 244.7 | 212.4 | 5 | 768 / 726 |
| re-cost arm (pixel radius 3, Q-max 8, C++ rule) | 329.4 | 251.9 | 77.5 | 1.31 | 324.6 | 250.5 | 5 | 768 / 726 |
| final (triangle, Q-max 8, C++ rule) | 322.4 | 247.0 | 75.4 | 1.31 | 319.1 | 246.1 | 5 | 735 / 693 |
<!-- /table:s1_per_call -->

<!-- table:s1_search -->
| arm | wall (s) | evaluate calls | us per evaluate | evaluate share (%) | cost | orientation sha1 | same in all rounds |
|---|---|---|---|---|---|---|---|
| base commit | 7.5 | 8311 | 800 | 88.5 | 0.4376801801801802 | 7a8743d0e2d1 | True |
| S1 | 5.5 | 8311 | 564 | 84.9 | 0.4376801801801802 | 7a8743d0e2d1 | True |
<!-- /table:s1_search -->

<!-- table:s1_stage_split -->
| part | before (us) | after (us) | saved (us) |
|---|---|---|---|
| evaluate() total (Q-max 8, 110 valid peaks) | 332.1 | 258.4 | 73.7 |
| stage A: loop over the peaks / frame_of_omegas | 19.0 | 1.7 | 17.3 |
| frame list: 174 get_image calls / table check | 30.9 | 1.0 | 29.9 |
| Bragg solve on g: original function / hoisted constants | 37.9 | 26.5 | 11.4 |
| stage D C call: compact list of 174 / table of 360 entries | 8.5 | 8.7 | -0.2 |
<!-- /table:s1_stage_split -->
