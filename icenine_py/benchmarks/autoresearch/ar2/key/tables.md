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
| voxel | source | err_deg | result_cost | cost_at_truth | carrier | carrier_err | carrier_rank | carrier_cost |
|---|---|---|---|---|---|---|---|---|
| 7194 | neighbor | 51.9 | 0.123 | 0.046 | 983 | 0.005 | 1 | 0.069 to 0.054 |
| 8890 | neighbor | 42.3 | 0.157 | 0.052 | 5876 | 0.023 | 1 | 0.123 to 0.078 |
| 14299 | neighbor | 44.5 | 0.137 | 0.053 | 9015 | 0.017 | 1 | 0.113 to 0.056 |
| 16508 | revisit | 33.5 | 0.192 | 0.057 | 8948 | 0.009 | 2 | 0.151 to 0.119 |
| 16534 | neighbor | 29.6 | 0.135 | 0.045 | 9015 | 0.017 | 1 | 0.113 to 0.056 |
| 17286 | neighbor | 31.4 | 0.121 | 0.024 | 6684 | 0.015 | 1 | 0.046 to 0.046 |
| 21180 | neighbor | 33.7 | 0.143 | 0.052 | 21170 | 0.004 | 1 | 0.086 to 0.066 |
<!-- /table:ar2b_break_trace -->

<!-- table:ar2b_break_trace_fused -->
| voxel | source | err_deg | result_cost | cost_at_truth | carrier | carrier_err | carrier_rank | carrier_cost |
|---|---|---|---|---|---|---|---|---|
| 1402 | neighbor | 46.0 | 0.111 | 0.051 | 19480 | 0.003 | 1 | 0.088 to 0.029 |
| 7179 | neighbor | 30.7 | 0.115 | 0.051 | 19480 | 0.003 | 1 | 0.088 to 0.029 |
| 7194 | neighbor | 51.9 | 0.123 | 0.046 | 983 | 0.005 | 1 | 0.069 to 0.054 |
| 8866 | revisit | 50.9 | 0.134 | 0.053 | 8889 | 0.012 | 1 | 0.063 to 0.063 |
| 8890 | neighbor | 42.3 | 0.157 | 0.052 | 5876 | 0.023 | 1 | 0.123 to 0.078 |
| 14299 | neighbor | 44.5 | 0.137 | 0.053 | 9015 | 0.017 | 1 | 0.113 to 0.056 |
| 16508 | revisit | 33.5 | 0.192 | 0.057 | 8948 | 0.009 | 2 | 0.151 to 0.119 |
| 16534 | neighbor | 29.6 | 0.135 | 0.045 | 9015 | 0.017 | 1 | 0.113 to 0.056 |
| 17286 | neighbor | 31.4 | 0.121 | 0.024 | 6684 | 0.015 | 1 | 0.046 to 0.046 |
| 21180 | neighbor | 33.7 | 0.143 | 0.052 | 21170 | 0.004 | 1 | 0.086 to 0.066 |
| 21234 | revisit | 46.8 | 0.118 | 0.057 | 3205 | 0.003 | 2 | 0.181 to 0.036 |
| 21328 | revisit | 47.2 | 0.131 | 0.045 | 9027 | 0.021 | 1 | 0.111 to 0.053 |
<!-- /table:ar2b_break_trace_fused -->

<!-- table:ar2b_breaks -->
| voxel | kind | source_on | err_on | taken | crossing | origin | against_path |
|---|---|---|---|---|---|---|---|
| 7194 | edge | neighbor | 51.86 | grain 154 | yes | carried from a right full search (crossing) | no |
| 8890 | edge | neighbor | 42.27 | grain 148 | yes | carried from a right full search (crossing) | no |
| 14299 | edge | neighbor | 44.52 | grain 24 | yes | carried from a right full search (crossing) | no |
| 16508 | edge | revisit | 33.49 | grain 400 | yes | carried from a right full search (crossing) | no |
| 16534 | edge | neighbor | 29.65 | grain 24 | yes | carried from a right full search (crossing) | no |
| 17286 | edge | neighbor | 31.44 | grain 247 | yes | carried from a right full search (crossing) | no |
| 21180 | edge | neighbor | 33.73 | grain 396 | yes | carried from a right full search (crossing) | no |
<!-- /table:ar2b_breaks -->

<!-- table:ar2b_breaks_fused -->
| voxel | kind | source_on | err_on | taken | crossing | origin | against_path |
|---|---|---|---|---|---|---|---|
| 1402 | edge | neighbor | 46.00 | grain 350 | yes | carried from a right full search (crossing) | no |
| 7179 | edge | neighbor | 30.74 | grain 350 | yes | carried from a right full search (crossing) | no |
| 7194 | edge | neighbor | 51.86 | grain 154 | yes | carried from a right full search (crossing) | no |
| 8866 | edge | revisit | 50.86 | grain 398 | yes | carried from a right full search (crossing) | no |
| 8890 | edge | neighbor | 42.27 | grain 148 | yes | carried from a right full search (crossing) | no |
| 14299 | edge | neighbor | 44.52 | grain 24 | yes | carried from a right full search (crossing) | no |
| 16508 | edge | revisit | 33.49 | grain 400 | yes | carried from a right full search (crossing) | no |
| 16534 | edge | neighbor | 29.65 | grain 24 | yes | carried from a right full search (crossing) | no |
| 17286 | edge | neighbor | 31.44 | grain 247 | yes | carried from a right full search (crossing) | no |
| 21180 | edge | neighbor | 33.73 | grain 396 | yes | carried from a right full search (crossing) | no |
| 21234 | edge | revisit | 46.81 | grain 300 | yes | carried from a right full search (crossing) | no |
| 21328 | edge | revisit | 47.23 | grain 151 | yes | carried from a right full search (crossing) | no |
<!-- /table:ar2b_breaks_fused -->

<!-- table:ar2b_candidates -->
no wrong indexed result behind a break
<!-- /table:ar2b_candidates -->

<!-- table:ar2b_evals -->
| run | total | full | local | full_searches | local_fits |
|---|---|---|---|---|---|
| key off (stored rerun) | 25,547,935 | 22,882,865 | 2,665,070 | 69 | 4049 |
| key on | 2,819,906 | 122,198 | 2,697,708 | 64 | 4095 |
| key off, FusedCost on | 25,602,383 | 22,904,126 | 2,698,257 | 69 | 4092 |
| key on, FusedCost on | 2,816,897 | 122,154 | 2,694,743 | 64 | 4095 |
<!-- /table:ar2b_evals -->

<!-- table:ar2b_flip_kinds -->
| flip | edge | corner | interior |
|---|---|---|---|
| fixed | 22 | 0 | 0 |
| broken | 7 | 0 | 0 |
<!-- /table:ar2b_flip_kinds -->

<!-- table:ar2b_index -->
| item | value |
|---|---|
| seeds and reseeds | 64 |
| indexed (every seed and reseed) | 64 |
| used (passed the seed test) | 64 |
| fallbacks to the grid search | 0 |
|   no candidate | 0 |
|   failed the seed test | 0 |
| evaluations of used indexing paths | 122,198 |
| evaluations of indexing paths that fell back | 0 |
| evaluations of the grid fallbacks | 0 |
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
| key on (FusedCost off) | 1272 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
| key on (FusedCost on) | 213 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
| key off (FusedCost on) | 1633 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
<!-- /table:ar2b_wall -->
