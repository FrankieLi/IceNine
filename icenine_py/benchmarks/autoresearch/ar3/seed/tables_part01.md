<!-- table:ar3f_check -->
| group | n | agree | disagree |
|---|---|---|---|
| fallback voxels (stored: fail) | 189 | 189 | 0 |
| control index seeds (stored: pass) | 200 | 200 | 0 |
<!-- /table:ar3f_check -->

<!-- table:ar3f_classes -->
| set | n | a | b | c | d | e | f |
|---|---|---|---|---|---|---|---|
| all fallback voxels | 189 | 0 | 1 | 0 | 175 | 13 | 0 |
| set A | 91 | 0 | 0 | 0 | 85 | 6 | 0 |
| set B | 98 | 0 | 1 | 0 | 90 | 7 | 0 |
<!-- /table:ar3f_classes -->

<!-- table:ar3f_cluster -->
| cluster | reseeds | fallbacks | evals_all | evals_fallbacks | evals_per_voxel |
|---|---|---|---|---|---|
| 1 | 260 | 36 | 12,900,978 | 12,315,955 | 49,619 |
| 2 | 101 | 11 | 3,984,483 | 3,755,994 | 19,725 |
| 3 to 5 | 130 | 10 | 3,724,448 | 3,414,981 | 7,570 |
| 6 or more | 463 | 17 | 7,055,934 | 5,808,028 | 118 |
<!-- /table:ar3f_cluster -->

<!-- table:ar3f_controls -->
| remedy | n | pass_s | stock_right | broken | fixed | add_median | add_mean | add_max |
|---|---|---|---|---|---|---|---|---|
| 1a more, best 2 | 200 | 200 of 200 (98 to 100%) | 193 | 0 | 0 | 15 | -1 | 281 |
| 1b more, best 5 | 200 | 200 of 200 (98 to 100%) | 193 | 1 | 3 | 1,912 | 2,042 | 4,447 |
| 2 cma-seed | 200 | 187 of 200 (89 to 96%) | 193 | 3 | 1 | -442 | -586 | 772 |
| 3 wide | 200 | 200 of 200 (98 to 100%) | 193 | 0 | 0 | 0 | 0 | 0 |
| 4 neighbours | 200 | 200 of 200 (98 to 100%) | 193 | 6 | 1 | -2,121 | -2,152 | 549 |
| 5 defer (tries only) | 200 | 187 of 200 (89 to 96%) | 193 | 6 | 1 | -2,137 | -2,312 | -869 |
| 6 outward | 200 | 200 of 200 (98 to 100%) | 193 | 0 | 0 | 0 | 0 | 0 |
| 7a order, first 2 (added) | 200 | 188 of 200 (90 to 97%) | 193 | 1 | 3 | -5 | 81 | 3,933 |
| 7b order, first 5 (added) | 200 | 200 of 200 (98 to 100%) | 193 | 1 | 4 | 2,024 | 2,162 | 6,224 |
<!-- /table:ar3f_controls -->

<!-- table:ar3f_gridend -->
| run | n | right_s | hit_ge_080 | median_err | with_nbr_within1 | median_nbr_within1 | median_nbr |
|---|---|---|---|---|---|---|---|
| D4 | 76 | 74 of 76 (91 to 99%) | 74 | 0.0046 | 76 | 7 | 27 |
| D | 142 | 141 of 142 (96 to 100%) | 140 | 0.0039 | 142 | 7 | 27 |
| P4 | 121 | 121 of 121 (97 to 100%) | 121 | 0.0032 | 121 | 7 | 27 |
| P | 172 | 172 of 172 (98 to 100%) | 172 | 0.0047 | 172 | 8 | 27 |
<!-- /table:ar3f_gridend -->

<!-- table:ar3f_levels -->
| level | voxels | at_boundary | inside |
|---|---|---|---|
| 0 | 1 | 1 | 0 |
| 0.1 | 1 | 1 | 0 |
| 0.2 | 3 | 3 | 0 |
| 0.5 | 34 | 34 | 0 |
| 1 | 150 | 150 | 0 |
| all | 189 | 189 | 0 |
| controls, level 0 | 14 | 11 | 3 |
| controls, level 0.1 | 9 | 7 | 2 |
| controls, level 0.2 | 21 | 18 | 3 |
| controls, level 0.5 | 36 | 31 | 5 |
| controls, level 1 | 120 | 111 | 9 |
| controls, all | 200 | 178 | 22 |
<!-- /table:ar3f_levels -->

<!-- table:ar3f_lists -->
| run | listed | counter | at_seed | at_reseed | at_paint |
|---|---|---|---|---|---|
| D4 | 76 | 76 | 2 | 74 | n/a |
| D | 142 | 142 | 4 | 138 | n/a |
| P4 | 121 | n/a | n/a | n/a | 121 |
| P | 172 | n/a | n/a | n/a | 172 |
<!-- /table:ar3f_lists -->

<!-- table:ar3f_pick -->
| chain | need_A | broken_A | fixed_A | total_A | need_B | broken_B | fixed_B | total_B |
|---|---|---|---|---|---|---|---|---|
| neighbours_only then more20_5 | 6 of 91 (3 to 14%) | 0 | 1 | 24,124 | 6 of 98 (3 to 13%) | 0 | 0 | 22,114 |
| neighbours then more20_5 | 6 of 91 (3 to 14%) | 0 | 1 | 24,779 | 6 of 98 (3 to 13%) | 0 | 0 | 22,553 |
| cma_seed then more20_5 | 6 of 91 (3 to 14%) | 0 | 1 | 28,587 | 10 of 98 (6 to 18%) | 0 | 0 | 41,001 |
| wide then more20_5 | 6 of 91 (3 to 14%) | 1 | 1 | 29,445 | 7 of 98 (4 to 14%) | 0 | 0 | 31,352 |
| more20_5 | 9 of 91 (5 to 18%) | 0 | 1 | 37,924 | 10 of 98 (6 to 18%) | 0 | 0 | 39,047 |
| more20_2 then more20_5 | 9 of 91 (5 to 18%) | 0 | 1 | 39,869 | 10 of 98 (6 to 18%) | 0 | 0 | 40,966 |
<!-- /table:ar3f_pick -->

<!-- table:ar3f_pick_extras -->
| chain | need_A | broken_A | fixed_A | total_A | need_B | broken_B | fixed_B | total_B |
|---|---|---|---|---|---|---|---|---|
| neighbours_only then order5 | 4 of 91 (2 to 11%) | 0 | 1 | 16,684 | 6 of 98 (3 to 13%) | 0 | 0 | 22,160 |
| neighbours then order5 | 4 of 91 (2 to 11%) | 0 | 1 | 17,339 | 6 of 98 (3 to 13%) | 0 | 0 | 22,599 |
| cma_seed then order5 | 4 of 91 (2 to 11%) | 0 | 1 | 21,199 | 9 of 98 (5 to 17%) | 0 | 0 | 37,679 |
| neighbours_only then order2 | 5 of 91 (2 to 12%) | 0 | 1 | 19,691 | 8 of 98 (4 to 15%) | 0 | 0 | 28,556 |
| neighbours then order2 | 5 of 91 (2 to 12%) | 0 | 1 | 20,347 | 8 of 98 (4 to 15%) | 0 | 0 | 28,995 |
| wide then order5 | 5 of 91 (2 to 12%) | 1 | 1 | 25,795 | 6 of 98 (3 to 13%) | 0 | 0 | 28,031 |
<!-- /table:ar3f_pick_extras -->

<!-- table:ar3f_rank -->
| set | n | none | rank0 | rank1_4 | rank5_19 | rank20_49 | median_cand | truth_hit_below | median_truth_hit | median_first_hit |
|---|---|---|---|---|---|---|---|---|---|---|
| fallback voxels | 189 | 1 | 145 | 43 | 0 | 0 | 50 | 0 | 1.000 | 0.435 |
| controls | 200 | 0 | 151 | 49 | 0 | 0 | 50 | 0 | 1.000 | 0.650 |
<!-- /table:ar3f_rank -->

<!-- table:ar3f_rem_A -->
| remedy | n | pass_s | need_grid | right_s | fixed | broken | both_right | both_wrong | ev_median | ev_mean | ev_max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1a more, best 2 | 91 | 0 of 91 (0 to 4%) | 91 | 90 of 91 (94 to 100%) | 0 | 0 | 90 | 1 | 1,844 | 1,944 | 3,795 |
| 1b more, best 5 | 91 | 82 of 91 (82 to 95%) | 9 | 91 of 91 (96 to 100%) | 1 | 0 | 90 | 0 | 4,257 | 4,307 | 5,985 |
| 2 cma-seed | 91 | 3 of 91 (1 to 9%) | 88 | 90 of 91 (94 to 100%) | 0 | 0 | 90 | 1 | 2,007 | 2,007 | 2,007 |
| 3 wide | 91 | 4 of 91 (2 to 11%) | 87 | 89 of 91 (92 to 99%) | 0 | 1 | 89 | 1 | 2,830 | 2,931 | 4,781 |
| 4 neighbours | 91 | 61 of 91 (57 to 76%) | 30 | 91 of 91 (96 to 100%) | 1 | 0 | 90 | 0 | 305 | 952 | 3,956 |
| 5 defer (tries only) | 91 | 61 of 91 (57 to 76%) | 30 | 91 of 91 (96 to 100%) | 1 | 0 | 90 | 0 | 305 | 297 | 610 |
| 6 outward | 91 | 0 of 91 (0 to 4%) | 91 | 90 of 91 (94 to 100%) | 0 | 0 | 90 | 1 | 2,493 | 2,594 | 4,444 |
| 7a order, first 2 (added) | 91 | 77 of 91 (76 to 91%) | 14 | 91 of 91 (96 to 100%) | 1 | 0 | 90 | 0 | 2,092 | 2,162 | 3,775 |
| 7b order, first 5 (added) | 91 | 84 of 91 (85 to 96%) | 7 | 91 of 91 (96 to 100%) | 1 | 0 | 90 | 0 | 4,307 | 4,395 | 6,161 |
<!-- /table:ar3f_rem_A -->

<!-- table:ar3f_rem_B -->
| remedy | n | pass_s | need_grid | right_s | fixed | broken | both_right | both_wrong | ev_median | ev_mean | ev_max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1a more, best 2 | 98 | 0 of 98 (0 to 4%) | 98 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 1,886 | 1,919 | 3,385 |
| 1b more, best 5 | 98 | 88 of 98 (82 to 94%) | 10 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 4,217 | 4,362 | 6,837 |
| 2 cma-seed | 98 | 1 of 98 (0 to 6%) | 97 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 2,007 | 2,007 | 2,007 |
| 3 wide | 98 | 4 of 98 (2 to 10%) | 94 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 2,872 | 2,903 | 4,371 |
| 4 neighbours | 98 | 76 of 98 (68 to 85%) | 22 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 305 | 745 | 3,736 |
| 5 defer (tries only) | 98 | 76 of 98 (68 to 85%) | 22 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 305 | 306 | 549 |
| 6 outward | 98 | 0 of 98 (0 to 4%) | 98 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 2,534 | 2,566 | 4,034 |
| 7a order, first 2 (added) | 98 | 81 of 98 (74 to 89%) | 17 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 2,139 | 2,215 | 4,163 |
| 7b order, first 5 (added) | 98 | 89 of 98 (83 to 95%) | 9 | 97 of 98 (94 to 100%) | 0 | 0 | 97 | 1 | 4,348 | 4,507 | 7,330 |
<!-- /table:ar3f_rem_B -->

<!-- table:ar3f_small -->
| group | n | ok_s | pass_right | hit_gap_median | stored_ev | tries_ev | saved | stored_ev_all |
|---|---|---|---|---|---|---|---|---|
| cluster 1, fallback | 36 | 1 of 36 (0 to 14%) | 24 of 36 (50 to 80%) | 0.126 | 345,168 | 305 | 344,863 | 12,315,955 |
| cluster 1, index seed used | 224 | 26 of 224 (8 to 16%) | 176 of 224 (73 to 83%) | 0.034 | 63,024 | 6,405 | 56,619 | 585,023 |
| cluster 2, fallback | 11 | 0 of 11 (0 to 26%) | 8 of 11 (43 to 90%) | 0.116 | 0 | 0 | 0 | 3,755,994 |
| cluster 2, index seed used | 90 | 8 of 90 (5 to 17%) | 69 of 90 (67 to 84%) | 0.032 | 21,452 | 2,440 | 19,012 | 228,489 |
<!-- /table:ar3f_small -->
