<!-- table:ar2b_17361 -->
| item | value |
|---|---|
| truth candidate: rank in the indexer's list | 1 |
| truth candidate: distance from the truth (deg) | 0.0449 |
| truth candidate: strict cost | 0.1623 |
| strict cost at the exact truth | 0.0949 |
| strict cost 0.02 deg from the truth | 0.1556 |
| strict cost 0.045 deg from the truth | 0.2367 |
| strict cost 0.1 deg from the truth | 0.4356 |
| finisher, total evaluations (truth candidate) | 924 |
| stage 1 (FindOptimal, MC blocks): evaluations | 97 |
| stage 1: blocks run, blocks accepted | 3, 0 |
| stage 1: first MC step and box (deg) | 0.132, 0.329 |
| stage 1: cost before, after | 0.1623, 0.1623 |
| stage 2 (VarianceMinimizing): evaluations | 826 |
| stage 2: cost before, after | 0.1623, 0.1623 |
| restart of the finisher from its end point: cost, moved (deg) | 0.1623, 0.0000 |
| ray Gauss-Newton from the end point: cost, accepted | 0.1283, True |
| ray Gauss-Newton: moved (deg), distance from the truth after (deg) | 0.0627, 0.0190 |
| ray Gauss-Newton: strict, smooth evaluations | 2, 14 |
<!-- /table:ar2b_17361 -->

<!-- table:ar2b_accuracy -->
| run | wrong | fixed | broken | both_wrong | unresolved |
|---|---|---|---|---|---|
| key off (stored rerun) | 37 of 2000 (1.3 to 2.5%) |  |  |  | 2 |
| key on | 22 of 2000 (0.7 to 1.7%) | 22 | 7 | 15 | 0 |
| key off, FusedCost on | 24 of 2000 (0.8 to 1.8%) |  |  |  | 1 |
| key on, FusedCost on | 22 of 2000 (0.7 to 1.7%) | 14 | 12 | 10 | 0 |
<!-- /table:ar2b_accuracy -->

<!-- table:ar2b_break_trace -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 7194 | neighbor | 51.9 | 0.123 | 0.046 | 983 | 0.005 | 1.7 | 1 | 0.069 to 0.054 |
| 8890 | neighbor | 42.3 | 0.157 | 0.052 | 5876 | 0.023 | 1.7 | 1 | 0.123 to 0.078 |
| 14299 | neighbor | 44.5 | 0.137 | 0.053 | 9015 | 0.017 | 1.5 | 1 | 0.113 to 0.056 |
| 16508 | revisit | 33.5 | 0.192 | 0.057 | 8948 | 0.009 | 2.5 | 2 | 0.151 to 0.119 |
| 16534 | neighbor | 29.6 | 0.135 | 0.045 | 9015 | 0.017 | 2.5 | 1 | 0.113 to 0.056 |
| 17286 | neighbor | 31.4 | 0.121 | 0.024 | 6684 | 0.015 | 3.1 | 1 | 0.046 to 0.046 |
| 21180 | neighbor | 33.7 | 0.143 | 0.052 | 21170 | 0.004 | 2.6 | 1 | 0.086 to 0.066 |
<!-- /table:ar2b_break_trace -->

<!-- table:ar2b_break_trace_fused -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 1402 | neighbor | 46.0 | 0.111 | 0.051 | 19480 | 0.003 | 4.7 | 1 | 0.088 to 0.029 |
| 7179 | neighbor | 30.7 | 0.115 | 0.051 | 19480 | 0.003 | 1.2 | 1 | 0.088 to 0.029 |
| 7194 | neighbor | 51.9 | 0.123 | 0.046 | 983 | 0.005 | 1.7 | 1 | 0.069 to 0.054 |
| 8866 | revisit | 50.9 | 0.134 | 0.053 | 8889 | 0.012 | 3.6 | 1 | 0.063 to 0.063 |
| 8890 | neighbor | 42.3 | 0.157 | 0.052 | 5876 | 0.023 | 1.7 | 1 | 0.123 to 0.078 |
| 14299 | neighbor | 44.5 | 0.137 | 0.053 | 9015 | 0.017 | 1.5 | 1 | 0.113 to 0.056 |
| 16508 | revisit | 33.5 | 0.192 | 0.057 | 8948 | 0.009 | 2.5 | 2 | 0.151 to 0.119 |
| 16534 | neighbor | 29.6 | 0.135 | 0.045 | 9015 | 0.017 | 2.5 | 1 | 0.113 to 0.056 |
| 17286 | neighbor | 31.4 | 0.121 | 0.024 | 6684 | 0.015 | 3.1 | 1 | 0.046 to 0.046 |
| 21180 | neighbor | 33.7 | 0.143 | 0.052 | 21170 | 0.004 | 2.6 | 1 | 0.086 to 0.066 |
| 21234 | revisit | 46.8 | 0.118 | 0.057 | 3205 | 0.003 | 7.0 | 2 | 0.181 to 0.036 |
| 21328 | revisit | 47.2 | 0.131 | 0.045 | 9027 | 0.021 | 1.2 | 1 | 0.111 to 0.053 |
<!-- /table:ar2b_break_trace_fused -->

<!-- table:ar2b_breaks -->
| voxel | kind | source_on | source_off | err_on | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 7194 | edge | neighbor | revisit | 51.86 | grain 154 | yes | matches a right full search | no | no | no |
| 8890 | edge | neighbor | neighbor | 42.27 | grain 148 | yes | matches a right full search | no | no | no |
| 14299 | edge | neighbor | revisit | 44.52 | grain 24 | yes | matches a right full search | no | no | no |
| 16508 | edge | revisit | revisit | 33.49 | grain 400 | yes | matches a right full search | no | no | no |
| 16534 | edge | neighbor | revisit | 29.65 | grain 24 | yes | matches a right full search | no | no | no |
| 17286 | edge | neighbor | revisit | 31.44 | grain 247 | yes | matches a right full search | no | no | no |
| 21180 | edge | neighbor | neighbor | 33.73 | grain 396 | yes | matches a right full search | no | no | no |
<!-- /table:ar2b_breaks -->

<!-- table:ar2b_breaks_fused -->
| voxel | kind | source_on | source_off | err_on | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 1402 | edge | neighbor | revisit | 46.00 | grain 350 | yes | matches a right full search | no | no | no |
| 7179 | edge | neighbor | revisit | 30.74 | grain 350 | yes | matches a right full search | no | no | no |
| 7194 | edge | neighbor | revisit | 51.86 | grain 154 | yes | matches a right full search | no | no | no |
| 8866 | edge | revisit | revisit | 50.86 | grain 398 | yes | matches a right full search | no | no | no |
| 8890 | edge | neighbor | revisit | 42.27 | grain 148 | yes | matches a right full search | no | no | no |
| 14299 | edge | neighbor | revisit | 44.52 | grain 24 | yes | matches a right full search | no | no | no |
| 16508 | edge | revisit | seed | 33.49 | grain 400 | yes | matches a right full search | no | accepted | no |
| 16534 | edge | neighbor | revisit | 29.65 | grain 24 | yes | matches a right full search | no | no | no |
| 17286 | edge | neighbor | revisit | 31.44 | grain 247 | yes | matches a right full search | no | no | no |
| 21180 | edge | neighbor | revisit | 33.73 | grain 396 | yes | matches a right full search | no | no | no |
| 21234 | edge | revisit | revisit | 46.81 | grain 300 | yes | matches a right full search | no | no | no |
| 21328 | edge | revisit | revisit | 47.23 | grain 151 | yes | matches a right full search | no | no | no |
<!-- /table:ar2b_breaks_fused -->

<!-- table:ar2b_evals -->
| run | total | full | local | full_searches | local_fits |
|---|---|---|---|---|---|
| key off (stored rerun) | 25,547,935 | 22,882,865 | 2,665,070 | 69 | 4049 |
| key on | 2,819,906 | 122,198 | 2,697,708 | 64 | 4095 |
| key off, FusedCost on | 25,602,383 | 22,904,126 | 2,698,257 | 69 | 4092 |
| key on, FusedCost on | 2,816,897 | 122,154 | 2,694,743 | 64 | 4095 |
<!-- /table:ar2b_evals -->

<!-- table:ar2b_fix_trace -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 1384 | neighbor | 20.5 | 0.132 | 0.071 | 5872 | 0.001 | 5.3 | n/a | n/a |
| 1398 | neighbor | 49.5 | 0.124 | 0.060 | 19440 | 0.001 | 2.5 | n/a | n/a |
| 1981 | neighbor | 34.2 | 0.133 | 0.074 | 19440 | 0.001 | 2.5 | n/a | n/a |
| 3204 | revisit | 34.2 | 0.125 | 0.074 | 19440 | 0.001 | 5.8 | n/a | n/a |
| 5008 | revisit | 31.5 | 0.110 | 0.043 | 21243 | 0.005 | 2.1 | n/a | n/a |
| 5876 | unresolved | 42.3 | 0.342 | 0.054 | 8866 | 50.850 | 2.6 | n/a | n/a |
| 8880 | neighbor | 33.7 | 0.124 | 0.066 | 8866 | 50.850 | 3.0 | n/a | n/a |
| 8899 | revisit | 34.2 | 0.120 | 0.074 | 19440 | 0.001 | 4.9 | n/a | n/a |
| 8909 | neighbor | 44.1 | 0.127 | 0.048 | 21210 | 0.015 | 2.3 | n/a | n/a |
| 8940 | neighbor | 31.5 | 0.109 | 0.043 | 21243 | 0.005 | 1.2 | n/a | n/a |
| 8978 | revisit | 57.3 | 0.139 | 0.057 | 14283 | 0.006 | 2.1 | n/a | n/a |
| 8988 | neighbor | 52.0 | 0.152 | 0.071 | 21289 | 0.004 | 0.6 | n/a | n/a |
| 13672 | neighbor | 36.7 | 0.134 | 0.060 | 5872 | 0.001 | 4.0 | n/a | n/a |
| 15749 | seed | 45.8 | 0.150 | 0.011 | itself | 45.753 | 0.0 | n/a | n/a |
| 15907 | neighbor | 46.0 | 0.137 | 0.049 | 5872 | 0.001 | 8.1 | n/a | n/a |
| 16470 | neighbor | 45.2 | 0.157 | 0.055 | 24281 | 0.018 | 1.0 | n/a | n/a |
| 18179 | neighbor | 29.8 | 0.113 | 0.045 | 8932 | 0.009 | 4.0 | n/a | n/a |
| 19418 | neighbor | 36.7 | 0.133 | 0.060 | 5872 | 0.001 | 3.1 | n/a | n/a |
| 19436 | neighbor | 53.3 | 0.128 | 0.051 | 777 | 0.006 | 3.2 | n/a | n/a |
| 21188 | neighbor | 47.1 | 0.134 | 0.074 | 8866 | 50.850 | 4.0 | n/a | n/a |
| 21215 | neighbor | 54.5 | 0.130 | 0.044 | 21210 | 0.015 | 1.0 | n/a | n/a |
| 24334 | unresolved | 50.8 | 0.423 | 0.048 | 7944 | 0.009 | 1.5 | n/a | n/a |
<!-- /table:ar2b_fix_trace -->

<!-- table:ar2b_fix_trace_fused -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 1384 | neighbor | 20.5 | 0.131 | 0.071 | 5872 | 0.003 | 5.3 | n/a | n/a |
| 1964 | neighbor | 43.8 | 0.124 | 0.028 | 21138 | 0.017 | 2.1 | n/a | n/a |
| 3213 | neighbor | 54.0 | 0.137 | 0.036 | 9015 | 0.022 | 4.6 | n/a | n/a |
| 8879 | revisit | 47.1 | 0.156 | 0.052 | 8899 | 0.004 | 5.3 | n/a | n/a |
| 8978 | revisit | 57.3 | 0.142 | 0.057 | 14283 | 0.012 | 2.1 | n/a | n/a |
| 8988 | neighbor | 52.0 | 0.152 | 0.071 | 21289 | 0.006 | 0.6 | n/a | n/a |
| 13250 | neighbor | 46.3 | 0.137 | 0.036 | 5876 | 0.006 | 2.5 | n/a | n/a |
| 13672 | neighbor | 36.7 | 0.136 | 0.060 | 5872 | 0.003 | 4.0 | n/a | n/a |
| 15907 | neighbor | 46.0 | 0.139 | 0.049 | 5872 | 0.003 | 8.1 | n/a | n/a |
| 16476 | unresolved | 49.1 | 0.387 | 0.028 | 399 | 0.012 | 3.1 | n/a | n/a |
| 16543 | neighbor | 22.1 | 0.115 | 0.086 | 16539 | 0.003 | 2.0 | n/a | n/a |
| 19418 | neighbor | 36.7 | 0.136 | 0.060 | 5872 | 0.003 | 3.1 | n/a | n/a |
| 19436 | neighbor | 53.3 | 0.131 | 0.051 | 777 | 0.013 | 3.2 | n/a | n/a |
| 21164 | revisit | 50.8 | 0.139 | 0.052 | 21153 | 0.029 | 3.0 | n/a | n/a |
<!-- /table:ar2b_fix_trace_fused -->

<!-- table:ar2b_fixes -->
| voxel | kind | source_off | source_on | err_off | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 1384 | edge | neighbor | revisit | 20.54 | grain 75 | yes | matches a right full search | no | no | no |
| 1398 | edge | neighbor | revisit | 49.49 | grain 23 | yes | matches a right full search | no | no | no |
| 1981 | edge | neighbor | revisit | 34.20 | grain 23 | yes | matches a right full search | no | no | no |
| 3204 | edge | revisit | revisit | 34.20 | grain 23 | yes | matches a right full search | no | no | no |
| 5008 | edge | revisit | neighbor | 31.45 | grain 399 | yes | matches a right full search | no | no | no |
| 5876 | edge | unresolved | seed | 42.25 | grain 398 | yes | matches a wrong grid seed | grid | rejected | indexed |
| 8880 | edge | neighbor | neighbor | 33.73 | grain 398 | yes | matches a wrong grid seed | grid | no | no |
| 8899 | edge | revisit | revisit | 34.20 | grain 23 | yes | matches a right full search | no | no | no |
| 8909 | edge | neighbor | neighbor | 44.11 | grain 332 | yes | matches a right full search | no | no | no |
| 8940 | edge | neighbor | neighbor | 31.45 | grain 399 | yes | matches a right full search | no | no | no |
| 8978 | edge | revisit | revisit | 57.34 | grain 401 | yes | matches a right full search | no | no | no |
| 8988 | edge | neighbor | revisit | 52.01 | grain 153 | yes | matches a right full search | no | no | no |
| 13672 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 15749 | edge | seed | revisit | 45.75 | grain 331 | yes | own grid seed (wrong) | grid | accepted | no |
| 15907 | edge | neighbor | revisit | 46.02 | grain 75 | yes | matches a right full search | no | no | no |
| 16470 | edge | neighbor | revisit | 45.18 | grain 145 | yes | matches a right full search | no | no | no |
| 18179 | edge | neighbor | neighbor | 29.78 | grain 300 | yes | matches a right full search | no | no | no |
| 19418 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 19436 | edge | neighbor | revisit | 53.33 | grain 279 | yes | matches a right full search | no | no | no |
| 21188 | edge | neighbor | revisit | 47.14 | grain 398 | yes | matches a wrong grid seed | grid | no | no |
| 21215 | edge | neighbor | neighbor | 54.51 | grain 332 | yes | matches a right full search | no | no | no |
| 24334 | edge | unresolved | seed | 50.85 | grain 56 | yes | matches a right full search | no | rejected | indexed |
<!-- /table:ar2b_fixes -->

<!-- table:ar2b_fixes_fused -->
| voxel | kind | source_off | source_on | err_off | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 1384 | edge | neighbor | revisit | 20.54 | grain 75 | yes | matches a right full search | no | no | no |
| 1964 | edge | neighbor | neighbor | 43.79 | grain 395 | yes | matches a right full search | no | no | no |
| 3213 | edge | neighbor | revisit | 54.04 | grain 24 | yes | matches a right full search | no | no | no |
| 8879 | edge | revisit | revisit | 47.13 | grain 147 | yes | matches a right full search | no | no | no |
| 8978 | edge | revisit | revisit | 57.34 | grain 401 | yes | matches a right full search | no | no | no |
| 8988 | edge | neighbor | revisit | 52.01 | grain 153 | yes | matches a right full search | no | no | no |
| 13250 | edge | neighbor | revisit | 46.25 | grain 148 | yes | matches a right full search | no | no | no |
| 13672 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 15907 | edge | neighbor | revisit | 46.02 | grain 75 | yes | matches a right full search | no | no | no |
| 16476 | edge | unresolved | seed | 49.08 | grain 144 | yes | matches a right full search | no | rejected | indexed |
| 16543 | edge | neighbor | revisit | 22.10 | grain 155 | yes | matches a right full search | no | no | no |
| 19418 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 19436 | edge | neighbor | revisit | 53.33 | grain 279 | yes | matches a right full search | no | no | no |
| 21164 | edge | revisit | revisit | 50.83 | grain 397 | yes | matches a right full search | no | no | no |
<!-- /table:ar2b_fixes_fused -->

<!-- table:ar2b_flip_kinds -->
| flip | edge | corner | interior |
|---|---|---|---|
| fixed | 22 | 0 | 0 |
| broken | 7 | 0 | 0 |
<!-- /table:ar2b_flip_kinds -->

<!-- table:ar2b_flips -->
| pair | flip | traced_in | n | crossing | right_seed_behind | wrong_grid_seed_behind | wrong_indexed_seed_behind | no_match | grid_seed_off | indexed_seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| stored pair | fixed | key off | 22 | 22 | 18 | 4 | 0 | 0 | 3 (2 rejected) | 2 |
| stored pair | broken | key on | 7 | 7 | 7 | 0 | 0 | 0 | 0 (0 rejected) | 0 |
| FusedCost-on pair | fixed | key off | 14 | 14 | 14 | 0 | 0 | 0 | 1 (1 rejected) | 1 |
| FusedCost-on pair | broken | key on | 12 | 12 | 12 | 0 | 0 | 0 | 1 (0 rejected) | 0 |
<!-- /table:ar2b_flips -->

<!-- table:ar2b_fused_diff -->
| pair | orientation_differs | median_deg | above_0_001_deg | above_0_01_deg | max_deg | mic_lines_differ | evals_differ |
|---|---|---|---|---|---|---|---|
| key on: FusedCost on against off | 2 of 2000 | 1.9e-02 | 2 | 2 | 0.019 | 45 of 2000 | 431 of 2000 |
| key off: FusedCost on against off | 2000 of 2000 | 1.3e-02 | 2000 | 1303 | 54.511 | 2000 of 2000 | 829 of 2000 |
<!-- /table:ar2b_fused_diff -->

<!-- table:ar2b_index -->
| item | value |
|---|---|
| seeds and reseeds | 64 |
| indexed (every seed and reseed) | 64 |
| used (passed the seed test) | 64 |
| fallbacks to the grid search | 0 of 64 (0.0 to 5.7%) |
|   no candidate | 0 |
|   failed the seed test | 0 |
| evaluations of used indexing paths | 122,198 |
| evaluations of indexing paths that fell back | 0 |
| evaluations of the grid fallbacks | 0 |
| key off: grid seeds and reseeds that failed the seed test | 7 of 69 (5.0 to 19.5%) |
| key on: indexed seeds and reseeds that failed the seed test | 0 of 64 (0.0 to 5.7%) |
| key off: voxels whose result is their own grid search, wrong | 2 of 62 (0.9 to 11.0%) |
| key on: voxels whose result is their own indexed seed, wrong | 0 of 64 (0.0 to 5.7%) |
<!-- /table:ar2b_index -->

<!-- table:ar2b_kinds -->
| run | edge | corner | interior |
|---|---|---|---|
| key off (stored rerun) | 37 | 0 | 0 |
| key on | 22 | 0 | 0 |
<!-- /table:ar2b_kinds -->

<!-- table:ar2b_noise -->
| pair | wrong_a | wrong_b | fixed | broken | both_wrong |
|---|---|---|---|---|---|
| key off with FusedCost on, against key off without | 37 | 24 | 20 | 7 | 17 |
<!-- /table:ar2b_noise -->

<!-- table:ar2b_times -->
| item | value |
|---|---|
| one-time pairing of the rays (s) | 17.6 |
| indexing, total (s) | 11.9 |
| indexing per voxel, mean (s) | 0.186 |
| indexing per voxel, median (s) | 0.181 |
| indexing per voxel, max (s) | 0.333 |
| whole path per voxel, mean (s) | 1.01 |
<!-- /table:ar2b_times -->

<!-- table:ar2b_wall -->
| run | wall_s | conditions |
|---|---|---|
| key off, stored rerun (FusedCost off) | 8997 | contended: 3 BFS runs at once, 1 thread each, 2026-10-29 |
| key on (FusedCost off) | 1272 | contended: timing preflight refused (an editor process at about 200% CPU); not a single run: overlapped with the voxel 17361 run (its log spans 23 s) and with short checks whose times were not logged |
| key on (FusedCost on) | 213 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
| key off (FusedCost on) | 1633 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
<!-- /table:ar2b_wall -->
