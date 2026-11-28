<!-- table:ar3b_accuracy -->
| region | n | C_p50 | C_p90 | C_p99 | D_p50 | D_p90 | D_p99 |
|---|---|---|---|---|---|---|---|
| layer | 24411 | 0.0172 | 0.0374 | 0.0566 | 0.0205 | 0.0446 | 0.0858 |
| pilot | 1980 | 0.0126 | 0.0349 | 0.0571 | 0.0142 | 0.0351 | 0.0814 |
| rest | 22431 | 0.0177 | 0.0374 | 0.0566 | 0.0207 | 0.0451 | 0.0858 |
<!-- /table:ar3b_accuracy -->

<!-- table:ar3b_direct -->
| region | n | lowest_cost | highest_hit | differ |
|---|---|---|---|---|
| layer | 363 | 350 of 363 (94.0 to 97.9%) | 360 of 363 (97.6 to 99.7%) | 14 of 363 (2.3 to 6.4%) |
| layer, edge | 363 | 350 of 363 (94.0 to 97.9%) |  |  |
| layer, corner | 0 | 0 of 0 |  |  |
| layer, interior | 0 | 0 of 0 |  |  |
| pilot | 47 | 46 of 47 (88.9 to 99.6%) | 47 of 47 (92.4 to 100.0%) | 1 of 47 (0.4 to 11.1%) |
| pilot, edge | 47 | 46 of 47 (88.9 to 99.6%) |  |  |
| pilot, corner | 0 | 0 of 0 |  |  |
| pilot, interior | 0 | 0 of 0 |  |  |
| rest | 316 | 304 of 316 (93.5 to 97.8%) | 313 of 316 (97.2 to 99.7%) | 13 of 316 (2.4 to 6.9%) |
| rest, edge | 316 | 304 of 316 (93.5 to 97.8%) |  |  |
| rest, corner | 0 | 0 of 0 |  |  |
| rest, interior | 0 | 0 of 0 |  |  |
<!-- /table:ar3b_direct -->

<!-- table:ar3b_evals -->
| run | total | seed | score | fit | grid | pilot | rest | wrong |
|---|---|---|---|---|---|---|---|---|
| C paint (key PaintGrains) | 1,316,345 | 990,350 | 79,979 | 246,016 | 0 | 124,324 | 1,192,021 | 14 of 24570 (0.03 to 0.10%) |
| D key-on BFS | 33,040,126 | 1,131,396 | 0 | 31,908,730 | 0 | 2,829,045 | 30,211,081 | 153 of 24570 (0.53 to 0.73%) |
<!-- /table:ar3b_evals -->

<!-- table:ar3b_fitlog_B -->
| moved | fits | share | evals | evals_share |
|---|---|---|---|---|
| not at all | 2,756 | 0.673 | 1,830,998 | 0.679 |
| under 0.01 | 932 | 0.228 | 619,494 | 0.230 |
| 0.01 to 0.1 | 6 | 0.001 | 3,752 | 0.001 |
| over 0.1 | 401 | 0.098 | 240,499 | 0.089 |
<!-- /table:ar3b_fitlog_B -->

<!-- table:ar3b_fitlog_D -->
| moved | fits | share | evals | evals_share |
|---|---|---|---|---|
| not at all | 27,644 | 0.568 | 18,345,781 | 0.575 |
| under 0.01 | 15,501 | 0.319 | 10,317,886 | 0.323 |
| 0.01 to 0.1 | 91 | 0.002 | 51,420 | 0.002 |
| over 0.1 | 5,403 | 0.111 | 3,193,643 | 0.100 |
<!-- /table:ar3b_fitlog_D -->

<!-- table:ar3b_fitlog_D_regions -->
| moved | pilot | rest | pilot_share | rest_share |
|---|---|---|---|---|
| not at all | 2,175 | 25,469 | 0.521 | 0.573 |
| under 0.01 | 1,533 | 13,968 | 0.367 | 0.314 |
| 0.01 to 0.1 | 13 | 78 | 0.003 | 0.002 |
| over 0.1 | 453 | 4,950 | 0.109 | 0.111 |
<!-- /table:ar3b_fitlog_D_regions -->

<!-- table:ar3b_fitlog_byvoxel -->
| run | region | voxels | moved_over_0p1 | none_over_0p1 | not_fitted |
|---|---|---|---|---|---|
| B | pilot | 2,000 | 4 of 353 (0.4 to 2.9%) | 18 of 1614 (0.7 to 1.8%) | 0 of 33 |
| D | layer | 24,570 | 25 of 4747 (0.36 to 0.78%) | 128 of 19491 (0.55 to 0.78%) | 0 of 332 |
| D | pilot | 2,000 | 6 of 401 (0.7 to 3.2%) | 13 of 1571 (0.5 to 1.4%) | 0 of 28 |
| D | rest | 22,570 | 19 of 4346 (0.28 to 0.68%) | 115 of 17920 (0.53 to 0.77%) | 0 of 304 |
<!-- /table:ar3b_fitlog_byvoxel -->

<!-- table:ar3b_kinds -->
| region | kind | voxels | C | D | crossC | crossD |
|---|---|---|---|---|---|---|
| layer | edge | 8,873 | 14 | 153 | 14 | 153 |
| layer | corner | 8,022 | 0 | 0 | 0 | 0 |
| layer | interior | 7,675 | 0 | 0 | 0 | 0 |
| pilot | edge | 813 | 1 | 19 | 1 | 19 |
| pilot | corner | 723 | 0 | 0 | 0 | 0 |
| pilot | interior | 464 | 0 | 0 | 0 | 0 |
| rest | edge | 8,060 | 13 | 134 | 13 | 134 |
| rest | corner | 7,299 | 0 | 0 | 0 | 0 |
| rest | interior | 7,211 | 0 | 0 | 0 | 0 |
<!-- /table:ar3b_kinds -->

<!-- table:ar3b_paired -->
| region | voxels | C | D | fixed | broken | both | intervals |
|---|---|---|---|---|---|---|---|
| layer | 24,570 | 14 of 24570 (0.03 to 0.10%) | 153 of 24570 (0.53 to 0.73%) | 145 | 6 | 8 | do not overlap |
| pilot | 2,000 | 1 of 2000 (0.01 to 0.28%) | 19 of 2000 (0.6 to 1.5%) | 19 | 1 | 0 | do not overlap |
| rest | 22,570 | 13 of 22570 (0.03 to 0.10%) | 134 of 22570 (0.50 to 0.70%) | 126 | 5 | 8 | do not overlap |
<!-- /table:ar3b_paired -->

<!-- table:ar3b_side -->
| m | voxels | evals | fixed | broken |
|---|---|---|---|---|
| 0.001 | 1 | 1,324 | 0 | 0 |
| 0.003 | 4 | 5,296 | 0 | 0 |
| 0.01 | 36 | 47,664 | 0 | 0 |
| 0.03 | 76 | 101,126 | 0 | 0 |
| 0.1 | 337 | 448,196 | 0 | 0 |
| 0.3 | 363 | 482,493 | 0 | 0 |
| 1 | 363 | 482,493 | 0 | 0 |
| all | 363 | 482,493 | 0 | 0 |
<!-- /table:ar3b_side -->
