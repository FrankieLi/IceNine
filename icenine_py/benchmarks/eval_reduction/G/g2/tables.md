<!-- table:evalG2_e0 -->
| n | wrong_base | wrong_integrated | fixed | broken | both_wrong | broken_traced | evals_base | evals_integrated | smooth | identical_to_task_e |
|---|---|---|---|---|---|---|---|---|---|---|
| 200 | 41 (15.5-26.6%) | 7 (1.7-7.0%) | 34 | 0 | 7 | 0 / 0 / 0 | 48700 | 49215 | 252 | 200 of 200 |
<!-- /table:evalG2_e0 -->

<!-- table:evalG2_ledger_stages -->
| stage | task_a | integrated | change |
|---|---|---|---|
| screen | 465057 | 466695 | +0.4% |
| recost | 402506 | 404295 | +0.4% |
| quick | 1098888 | 1110036 | +1.0% |
| find | 17312 | 15792 | -8.8% |
| variance | 8509 | 10093 | +18.6% |
| final | 6 | 6 | +0.0% |
| all | 1992278 | 2006917 | +0.7% |
| smooth (E3+E4, outside the strict count) | 0 | 90183 | n/a |
<!-- /table:evalG2_ledger_stages -->

<!-- table:evalG2_ledger_voxels -->
| voxel | evals_task_a | evals_key_off | evals_integrated | change | smooth | err_task_a | err_integrated | wall_key_off | wall_integrated |
|---|---|---|---|---|---|---|---|---|---|
| 15901 | 327619 | 327619 | 329464 | +0.6% | 14804 | 0.0378 | 0.0095 | 118 | 166 |
| 17034 | 340422 | 340422 | 342111 | +0.5% | 15231 | 0.0217 | 0.0210 | 127 | 170 |
| 18558 | 329915 | 329915 | 330488 | +0.2% | 15024 | 0.0050 | 0.0042 | 123 | 159 |
| 19889 | 329724 | 329724 | 334829 | +1.5% | 15232 | 0.0145 | 0.0278 | 118 | 180 |
| 22867 | 333046 | 333046 | 335507 | +0.7% | 14884 | 0.0119 | 0.0227 | 125 | 164 |
| 5242 | 331552 | 331552 | 334518 | +0.9% | 15008 | 0.0604 | 0.0116 | 118 | 162 |
<!-- /table:evalG2_ledger_voxels -->

<!-- table:evalG2_pilot150_accuracy -->
| run | n | wrong | unresolved | grains_lost | median_err_right | seeds | reseeds | neighbor_fits | accepted_first | revisits |
|---|---|---|---|---|---|---|---|---|---|---|
| keyoff_branch | 150 | 1 (0.1-3.7%) | 0 | 0 | 0.0139 | 4 | 4 (4 accepted) | 146 | 70 | 123 (72 accepted) |
| outward_only | 150 | 2 (0.4-4.7%) | 0 | 0 | 0.0209 | 4 | 5 (5 accepted) | 146 | 70 | 123 (71 accepted) |
| integrated | 150 | 2 (0.4-4.7%) | 0 | 0 | 0.0234 | 4 | 5 (5 accepted) | 146 | 70 | 123 (71 accepted) |
<!-- /table:evalG2_pilot150_accuracy -->

<!-- table:evalG2_pilot150_evals -->
| run | seed | reseed | neighbor | revisit | total | smooth_full_searches |
|---|---|---|---|---|---|---|
| keyoff_branch | 1318606 | 1327503 | 96026 | 81646 | 2823781 | 0 |
| outward_only | 1324683 | 1329474 | 1222221 | 790378 | 4666756 | 0 |
| integrated | 1338437 | 1344358 | 1264407 | 790356 | 4737558 | 0 |
<!-- /table:evalG2_pilot150_evals -->

<!-- table:evalG2_pilot150_paired -->
| a | b | wrong_a | wrong_b | fixed | broken | both_wrong |
|---|---|---|---|---|---|---|
| keyoff_branch | outward_only | 1 | 2 | 0 | 1 | 1 |
| keyoff_branch | integrated | 1 | 2 | 0 | 1 | 1 |
<!-- /table:evalG2_pilot150_paired -->

<!-- table:evalG2_pilot150_policy -->
| run | phase | runs | accepted | fallback | fallback_accepted | strict | screen | smooth | jac | strict_equiv | wall_s | median_wall_s | max_wall_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| outward_only | neighbor | 146 | 119 | 27 | 0 | 1453 | 1124739 | 20042 | 989 | 150717 | 99 | 0.13 | 17.7 |
| outward_only | reseed | 5 | 5 | 0 | 0 | 10 | 0 | 84 | 26 | 921 | 1 | 0.10 | 0.1 |
| outward_only | revisit | 123 | 106 | 17 | 0 | 977 | 708169 | 13417 | 768 | 101169 | 55 | 0.12 | 2.8 |
| integrated | neighbor | 146 | 119 | 28 | 1 | 1496 | 1166396 | 21074 | 1094 | 158919 | 85 | 0.15 | 2.8 |
| integrated | reseed | 5 | 5 | 0 | 0 | 10 | 0 | 98 | 31 | 1089 | 1 | 0.10 | 0.2 |
| integrated | revisit | 123 | 106 | 17 | 0 | 977 | 708169 | 13602 | 839 | 103507 | 56 | 0.12 | 2.8 |
<!-- /table:evalG2_pilot150_policy -->

<!-- table:evalG2_pilot150_wall -->
| run | seed | reseed | neighbor | revisit | total |
|---|---|---|---|---|---|
| keyoff_branch | 466 | 468 | 43 | 36 | 1013 |
| outward_only | 472 | 467 | 142 | 91 | 1171 |
| integrated | 653 | 638 | 129 | 92 | 1512 |
<!-- /table:evalG2_pilot150_wall -->
