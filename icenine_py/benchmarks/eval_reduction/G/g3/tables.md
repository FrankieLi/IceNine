<!-- table:evalG3_bfs_accuracy -->
| run | n | wrong | unresolved | grains_lost | median_err_right | seeds | reseeds | neighbor_fits | accepted_first | revisits |
|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 2000 | 37 (1.3-2.5%) | 1 | 0 | 0.0118 | 34 | 35 (32 accepted) | 1966 | 770 | 2126 (1165 accepted) |
| keyoff_rerun | 2000 | 37 (1.3-2.5%) | 2 | 0 | 0.0101 | 40 | 29 (28 accepted) | 1960 | 769 | 2089 (1167 accepted) |
| smooth_only | 2000 | 34 (1.2-2.4%) | 0 | 0 | 0.0082 | 35 | 30 (30 accepted) | 1965 | 768 | 2108 (1170 accepted) |
| outward_only | 2000 | 36 (1.3-2.5%) | 1 | 1 | 0.0113 | 34 | 36 (35 accepted) | 1966 | 764 | 2166 (1168 accepted) |
| integrated | 2000 | 42 (1.6-2.8%) | 0 | 1 | 0.0077 | 38 | 32 (32 accepted) | 1962 | 762 | 2140 (1173 accepted) |
<!-- /table:evalG3_bfs_accuracy -->

<!-- table:evalG3_bfs_evals -->
| run | seed | reseed | neighbor | revisit | total | smooth_full_searches |
|---|---|---|---|---|---|---|
| mc_clean | 11274936 | 11598290 | 1290496 | 1398836 | 25562558 | 0 |
| keyoff_rerun | 13264873 | 9617992 | 1293861 | 1371209 | 25547935 | 0 |
| smooth_only | 11743593 | 10066910 | 1289656 | 1390172 | 24490331 | 980594 |
| outward_only | 11265144 | 10415877 | 19987665 | 19950797 | 61619483 | 0 |
| integrated | 12747993 | 9493718 | 19736111 | 18843288 | 60821110 | 996303 |
<!-- /table:evalG3_bfs_evals -->

<!-- table:evalG3_bfs_paired -->
| a | b | wrong_a | wrong_b | fixed | broken | both_wrong |
|---|---|---|---|---|---|---|
| mc_clean | keyoff_rerun | 37 | 37 | 14 | 14 | 23 |
| mc_clean | smooth_only | 37 | 34 | 13 | 10 | 24 |
| mc_clean | outward_only | 37 | 36 | 9 | 8 | 28 |
| mc_clean | integrated | 37 | 42 | 11 | 16 | 26 |
| keyoff_rerun | smooth_only | 37 | 34 | 14 | 11 | 23 |
| keyoff_rerun | outward_only | 37 | 36 | 18 | 17 | 19 |
| keyoff_rerun | integrated | 37 | 42 | 15 | 20 | 22 |
| smooth_only | integrated | 34 | 42 | 5 | 13 | 29 |
<!-- /table:evalG3_bfs_paired -->

<!-- table:evalG3_bfs_policy -->
| run | phase | runs | accepted | fallback | fallback_accepted | strict | screen | smooth | jac | strict_equiv | wall_s | median_wall_s | max_wall_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| outward_only | neighbor | 1966 | 1521 | 460 | 15 | 23290 | 18674835 | 325889 | 12560 | 2385335 | 1239 | 0.10 | 17.8 |
| outward_only | reseed | 36 | 33 | 3 | 0 | 201 | 124971 | 2599 | 209 | 20436 | 10 | 0.11 | 2.5 |
| outward_only | revisit | 2166 | 1722 | 449 | 5 | 23444 | 18495783 | 327026 | 14220 | 2420169 | 1237 | 0.11 | 3.0 |
| integrated | neighbor | 1962 | 1523 | 455 | 16 | 23028 | 18424908 | 325119 | 13461 | 2391401 | 1218 | 0.11 | 3.0 |
| integrated | reseed | 32 | 30 | 2 | 0 | 150 | 83314 | 1836 | 163 | 14561 | 7 | 0.09 | 2.6 |
| integrated | revisit | 2140 | 1722 | 425 | 7 | 22283 | 17412731 | 312038 | 14757 | 2325773 | 1175 | 0.11 | 3.1 |
<!-- /table:evalG3_bfs_policy -->

<!-- table:evalG3_bfs_wall -->
| run | seed | reseed | neighbor | revisit | total |
|---|---|---|---|---|---|
| mc_clean | 6025 | 6213 | 723 | 783 | 13744 |
| keyoff_rerun | 4553 | 3295 | 559 | 590 | 8997 |
| smooth_only | 5381 | 4629 | 540 | 585 | 11136 |
| outward_only | 3852 | 3489 | 1791 | 1848 | 10980 |
| integrated | 5951 | 4325 | 1770 | 1775 | 13822 |
<!-- /table:evalG3_bfs_wall -->

<!-- table:evalG3_full_sample_estimate -->
| run | full_searches | full_search_wall_s | local_fit_wall_s | grain_scale | voxel_scale | estimate_h | timing |
|---|---|---|---|---|---|---|---|
| stored configuration, today's code | 69 | 7,848 | 1,149 | 8.28 | 12.29 | 22.0 | contended: 3 concurrent 2000-voxel BFS runs (Task G2) |
| SmoothCostRank only | 65 | 10,011 | 1,125 | 8.28 | 12.29 | 26.9 | contended: one 2000-voxel BFS run (Task G3); preflight quiet at the start, other jobs on the machine not controlled |
<!-- /table:evalG3_full_sample_estimate -->

<!-- table:evalG3_full_search_time -->
| run | full_searches | wall_per_full_search_s | strict_per_full_search | smooth_per_full_search | implied_ms_per_smooth | timing |
|---|---|---|---|---|---|---|
| stored configuration, today's code | 69 | 114 | 331,636 | 0 | n/a | contended: 3 concurrent 2000-voxel BFS runs (Task G2) |
| SmoothCostRank only | 65 | 154 | 335,546 | 15,086 | 2.58 | contended: one 2000-voxel BFS run (Task G3); preflight quiet at the start, other jobs on the machine not controlled |
| both (Task G2's integrated configuration) | 70 | 147 | 317,739 | 14,233 | 2.66 | contended: 2 concurrent 2000-voxel BFS runs plus the E0 run (2 workers) and a 150-voxel pilot at the start (Task G2) |
<!-- /table:evalG3_full_search_time -->

<!-- table:evalG3_levers -->
| lever | task | evaluations | wrong | time | carried |
|---|---|---|---|---|---|
| C++ overlap rule (`CppOverlapRule 1`) | A | 332,046 -> 228,146 per voxel (six full-sample voxels) | E0: 41 -> 56 (15 fixed / 30 broken) | 173 -> 116 s per voxel (mixed labels) | no (lenient rule stays, D-20261013-I6) |
| C++ edge rules, omega interpolation, batched screen (`CppEdgeRules`, `CppOmegaInterp`, `BatchedScreen`) | A | 332,046 -> 329,595 per voxel (six voxels, lenient rule, all three) | E0 (lenient rule, edge rules and interpolation): 41 -> 40 (21 fixed / 20 broken) | 122.5 -> 96.2 s, voxel 15901, C++ rule with `CppEdgeRules`, adding the other two (single-worker) | `BatchedScreen` on by default since 2026-10-13; the other two opt-in |
| candidate de-duplication, best arm (`CandidateDedup clique=cpp:1`) | C | 75.8% of the base (six voxels); 99.4% on E0 | E0: 41 -> 55 (10 fixed / 24 broken) | not measured single-worker | no (D-20261014-N4) |
| smooth ranking key of levels 0-2 (`SmoothCostRank E3+E4,tau=0.25`) | E, G1 | E0 48,700 -> 49,215 strict + 252 smooth per voxel; six voxels 332,046 -> 334,486 + 15,030 | E0: 41 -> 7 (34 fixed / 0 broken) | +21.1% per full-sample voxel with the cache (voxel 17034, single-worker) | yes (D-20261029-N1) |
| smooth ranking key of levels 0-2 (`SmoothCostRank E1+E2+E3+E4`) | E, G1 | E0 48,700 -> 49,653 strict + 263 smooth per voxel; six voxels 332,046 -> 348,018 + 15,755 | E0: 41 -> 5 (36 fixed / 0 broken) | +106.9% per full-sample voxel with the cache (voxel 17034, single-worker) | no (G1 default rule: fewer net fixes on the full sample) |
| SmoothCost cache per image set | G1 | none (the same results) | no change (bit-identical) | saves the 16.7 s build per voxel (17.6 / 17.2 s measured, single-worker) | yes, under the key |
| outward search for BFS neighbours and reseeds (`OutwardSearch 1`) | F, G2 | BFS region 25.5M -> 61.6M (mostly screen evaluations) | BFS region: 37 -> 36 (18 fixed / 17 broken) | 8,997 -> 10,980 s (contended) | no (D-20261029-N2) |
<!-- /table:evalG3_levers -->

<!-- table:evalG3_total_bfs -->
| run | wrong | fixed_broken | grains_lost | seed | reseed | neighbor | revisit | total | smooth | strict_equivalent | equiv_change | wall_s | timing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| stored configuration, today's code | 37 (1.3-2.5%) | - | 0 | 13,264,873 | 9,617,992 | 1,293,861 | 1,371,209 | 25,547,935 | 0 | 25,547,935 | +0.0% | 8,997 | contended: 3 concurrent 2000-voxel BFS runs (Task G2) |
| SmoothCostRank only | 34 (1.2-2.4%) | 14 / 11 | 0 | 11,743,593 | 10,066,910 | 1,289,656 | 1,390,172 | 24,490,331 | 980,594 | 28,703,943 | +12.4% | 11,136 | contended: one 2000-voxel BFS run (Task G3); preflight quiet at the start, other jobs on the machine not controlled |
| OutwardSearch only | 36 (1.3-2.5%) | 18 / 17 | 1 | 11,265,144 | 10,415,877 | 19,987,665 | 19,950,797 | 61,619,483 | 0 | 29,102,899 | +13.9% | 10,980 | contended: 2 concurrent 2000-voxel BFS runs plus the E0 run (2 workers) and a 150-voxel pilot at the start (Task G2) |
| both (Task G2's integrated configuration) | 42 (1.6-2.8%) | 15 / 20 | 1 | 12,747,993 | 9,493,718 | 19,736,111 | 18,843,288 | 60,821,110 | 996,303 | 33,867,545 | +32.6% | 13,822 | contended: 2 concurrent 2000-voxel BFS runs plus the E0 run (2 workers) and a 150-voxel pilot at the start (Task G2) |
<!-- /table:evalG3_total_bfs -->

<!-- table:evalG3_total_nostart_stages -->
| stage | before | after | change |
|---|---|---|---|
| screen | 465,057 | 466,695 | +0.4% |
| recost | 402,506 | 404,295 | +0.4% |
| quick | 1,098,888 | 1,110,036 | +1.0% |
| find | 17,312 | 15,792 | -8.8% |
| variance | 8,509 | 10,093 | +18.6% |
| final | 6 | 6 | +0.0% |
| strict, all stages | 1,992,278 | 2,006,917 | +0.7% |
| smooth key (E3+E4) | 0 | 90,183 | n/a |
| strict-equivalent (smooth x 4.297) | 1,992,278 | 2,394,433 | +20.2% |
<!-- /table:evalG3_total_nostart_stages -->

<!-- table:evalG3_total_sets -->
| set | n | wrong_before | wrong_after | fixed_broken | strict_before | strict_after | smooth | equiv_change |
|---|---|---|---|---|---|---|---|---|
| six full-sample voxels (sum) | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 / 0 | 1,992,278 | 2,006,917 | 90,183 | +20.2% |
| Task E's other 24 full-sample voxels (per voxel) | 24 | 1 (0.7-20.2%) | 0 (0.0-13.8%) | 1 / 0 | 330,584 | 335,047 | 15,119 | +21.0% |
| E0 200 voxels (per voxel) | 200 | 41 (15.5-26.6%) | 7 (1.7-7.0%) | 34 / 0 | 48,700 | 49,215 | 252 | +1.7% |
| BFS pilot region, 2,000 voxels (whole run) | 2000 | 37 (1.3-2.5%) | 34 (1.2-2.4%) | 14 / 11 | 25,547,935 | 24,490,331 | 980,594 | +12.4% |
<!-- /table:evalG3_total_sets -->

<!-- table:evalG3_trace -->
| a | b | flip | n | another grain | own full search | carried from a wrong full search | boundary crossing | untraced |
|---|---|---|---|---|---|---|---|---|
| stored configuration, today's code | SmoothCostRank only | broken (wrong in b) | 11 | 11 | 0 | 0 | 11 | 0 |
| stored configuration, today's code | SmoothCostRank only | fixed (wrong in a) | 14 | 14 | 1 | 3 | 10 | 0 |
| stored configuration, today's code | OutwardSearch only | broken (wrong in b) | 17 | 17 | 5 | 3 | 9 | 0 |
| stored configuration, today's code | OutwardSearch only | fixed (wrong in a) | 18 | 18 | 1 | 3 | 14 | 0 |
| stored configuration, today's code | both (Task G2's integrated configuration) | broken (wrong in b) | 20 | 20 | 4 | 4 | 12 | 0 |
| stored configuration, today's code | both (Task G2's integrated configuration) | fixed (wrong in a) | 15 | 15 | 1 | 3 | 11 | 0 |
| stored pilot (2026-10-10) | stored configuration, today's code | broken (wrong in b) | 14 | 14 | 1 | 3 | 10 | 0 |
| stored pilot (2026-10-10) | stored configuration, today's code | fixed (wrong in a) | 14 | 14 | 1 | 5 | 8 | 0 |
| SmoothCostRank only | both (Task G2's integrated configuration) | broken (wrong in b) | 13 | 13 | 2 | 4 | 7 | 0 |
| SmoothCostRank only | both (Task G2's integrated configuration) | fixed (wrong in a) | 5 | 5 | 0 | 0 | 5 | 0 |
<!-- /table:evalG3_trace -->

<!-- table:evalG3_trace_n13 -->
| flip | n | traced_to_arm | traced_to_draw | untraced |
|---|---|---|---|---|
| broken | 11 | 0 | 11 | 0 |
| fixed | 14 | 0 | 14 | 0 |
<!-- /table:evalG3_trace_n13 -->

<!-- table:evalG3_trace_probe -->
| full_sample_voxel | err_key_off | err_key_on | outcome |
|---|---|---|---|
| 8866 | 0.001 | 0.025 | both right |
| 15749 | 0.020 | 0.007 | both right |
<!-- /table:evalG3_trace_probe -->

<!-- table:evalG3_wrong_by_origin -->
| run | n | another grain | own full search | carried from a wrong full search | boundary crossing | untraced | max_distance_um | full_searches |
|---|---|---|---|---|---|---|---|---|
| stored pilot (2026-10-10) | 37 | 37 | 1 | 5 | 31 | 0 | 9.4 | 67 |
| stored configuration, today's code | 37 | 37 | 2 | 3 | 32 | 0 | 9.4 | 63 |
| SmoothCostRank only | 34 | 34 | 0 | 0 | 34 | 0 | 0.0 | 62 |
| OutwardSearch only | 36 | 36 | 5 | 4 | 27 | 0 | 9.4 | 68 |
| both (Task G2's integrated configuration) | 42 | 42 | 4 | 9 | 29 | 0 | 0.0 | 65 |
<!-- /table:evalG3_wrong_by_origin -->
