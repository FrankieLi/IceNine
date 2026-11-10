<!-- table:ar3c_control -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| this run P (paint) | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| this run P plus R | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| this run D (BFS) | 4,877 | 0.0169 | 0.0554 | 17 of 4877 (0.22 to 0.56%) | 17 of 4877 (0.22 to 0.56%) |
| step 2 C (paint) | 4,877 | 0.0177 | 0.0486 | 1 of 4877 (0.00 to 0.12%) | 1 of 4877 (0.00 to 0.12%) |
| step 2 D (BFS) | 4,877 | 0.0194 | 0.0486 | 37 of 4877 (0.6 to 1.0%) | 37 of 4877 (0.6 to 1.0%) |
<!-- /table:ar3c_control -->

<!-- table:ar3c_converged -->
| level | n | med_evals | share | med_share | start_within |
|---|---|---|---|---|---|
| level 0 | 4,774 | 662 | 0.999 | 0.998 | 4774 of 4774 (99.9 to 100.0%) |
| level 0.1 | 4,152 | 662 | 0.999 | 0.998 | 4152 of 4152 (99.9 to 100.0%) |
| level 0.2 | 4,332 | 662 | 0.998 | 0.998 | 4325 of 4332 (99.7 to 99.9%) |
| level 0.5 | 3,196 | 662 | 0.998 | 0.998 | 3186 of 3196 (99.4 to 99.8%) |
| level 1 | 2,559 | 662 | 0.998 | 0.998 | 2540 of 2559 (98.8 to 99.5%) |
| all | 19,013 | 662 | 0.998 | 0.998 | 18977 of 19013 (99.7 to 99.9%) |
<!-- /table:ar3c_converged -->

<!-- table:ar3c_dfits -->
| level | fits | evals | med_evals | mv_med | mv_p95 | same | s_med | e_med | e_p95 | e_wrong |
|---|---|---|---|---|---|---|---|---|---|---|
| level 0 | 8,617 | 5,695,834 | 662 | 0.0000 | 0.3226 | 4908 of 8617 (55.9 to 58.0%) | 0.0484 | 0.0484 | 53.2871 | 3846 of 8617 (43.6 to 45.7%) |
| level 0.1 | 11,893 | 7,875,695 | 662 | 0.0000 | 0.3705 | 6865 of 11893 (56.8 to 58.6%) | 0.1953 | 0.1954 | 54.0878 | 5930 of 11893 (49.0 to 50.8%) |
| level 0.2 | 18,608 | 12,369,093 | 662 | 0.0000 | 0.3983 | 10077 of 18608 (53.4 to 54.9%) | 20.7688 | 20.7151 | 54.3637 | 9843 of 18608 (52.2 to 53.6%) |
| level 0.5 | 22,516 | 15,059,166 | 662 | 0.0000 | 0.4481 | 11119 of 22516 (48.7 to 50.0%) | 0.4467 | 0.3998 | 53.4110 | 10847 of 22516 (47.5 to 48.8%) |
| level 1 | 28,558 | 19,320,043 | 662 | 0.0000 | 0.6240 | 12484 of 28558 (43.1 to 44.3%) | 0.6528 | 0.5059 | 53.8450 | 12963 of 28558 (44.8 to 46.0%) |
| all | 90,192 | 60,319,831 | 662 | 0.0000 | 0.5093 | 45453 of 90192 (50.1 to 50.7%) | 0.6252 | 0.4857 | 53.9183 | 43429 of 90192 (47.8 to 48.5%) |
<!-- /table:ar3c_dfits -->

<!-- table:ar3c_entries -->
| level | grains | entries | mean | mx | multi | none | grid |
|---|---|---|---|---|---|---|---|
| level 0 | 98 | 103 | 1.05 | 2 | 5 of 98 (2.2 to 11.4%) | 0 | 1 |
| level 0.1 | 81 | 334 | 4.12 | 7 | 81 of 81 (95.5 to 100.0%) | 0 | 1 |
| level 0.2 | 101 | 707 | 7.00 | 13 | 101 of 101 (96.3 to 100.0%) | 0 | 2 |
| level 0.5 | 101 | 1381 | 13.67 | 29 | 101 of 101 (96.3 to 100.0%) | 0 | 31 |
| level 1 | 116 | 2844 | 24.52 | 48 | 116 of 116 (96.8 to 100.0%) | 0 | 137 |
<!-- /table:ar3c_entries -->

<!-- table:ar3c_evals -->
| arm | evals |
|---|---|
| P (paint) | 78,142,641 |
| R (one fit on each painted voxel) | 12,843,129 |
| P plus R | 90,985,770 |
| R2 (second fit, every fourth voxel) | 3,217,669 |
| D (BFS, indexed seeds) | 118,350,225 |
<!-- /table:ar3c_evals -->

<!-- table:ar3c_evals_stage -->
| stage | evals |
|---|---|
| seed_index | 14,106,247 |
| score | 291,783 |
| fit_rejected | 3,459,521 |
| fit_accepted | 124,976 |
| index_rejected | 336,721 |
| seed_grid | 58,461,276 |
| grid_rejected | 1,362,117 |
| total | 78,142,641 |
<!-- /table:ar3c_evals_stage -->

<!-- table:ar3c_map_four -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| P plus R, all | 24,570 | 0.0342 | 0.1302 | 2531 of 24570 (9.9 to 10.7%) | 170 of 24570 (0.60 to 0.80%) |
| D, all | 24,570 | 0.0427 | 0.1459 | 3545 of 24570 (14.0 to 14.9%) | 203 of 24570 (0.72 to 0.95%) |
| P (start), all | 24,570 | 0.0342 | 0.1304 | 2541 of 24570 (10.0 to 10.7%) | 170 of 24570 (0.60 to 0.80%) |
| P plus R, level 0 | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| D, level 0 | 4,877 | 0.0169 | 0.0554 | 17 of 4877 (0.22 to 0.56%) | 17 of 4877 (0.22 to 0.56%) |
| P plus R, level 0.1 | 4,488 | 0.0313 | 0.0699 | 41 of 4488 (0.7 to 1.2%) | 10 of 4488 (0.12 to 0.41%) |
| D, level 0.1 | 4,488 | 0.0436 | 0.0964 | 179 of 4488 (3.5 to 4.6%) | 23 of 4488 (0.34 to 0.77%) |
| P plus R, level 0.2 | 5,041 | 0.0401 | 0.0897 | 152 of 5041 (2.6 to 3.5%) | 23 of 5041 (0.30 to 0.68%) |
| D, level 0.2 | 5,041 | 0.0533 | 0.1167 | 554 of 5041 (10.2 to 11.9%) | 38 of 5041 (0.5 to 1.0%) |
| P plus R, level 0.5 | 4,625 | 0.0570 | 0.1374 | 770 of 4625 (15.6 to 17.7%) | 39 of 4625 (0.6 to 1.2%) |
| D, level 0.5 | 4,625 | 0.0665 | 0.1618 | 1251 of 4625 (25.8 to 28.3%) | 47 of 4625 (0.8 to 1.3%) |
| P plus R, level 1 | 5,539 | 0.0550 | 0.1873 | 1565 of 5539 (27.1 to 29.5%) | 95 of 5539 (1.4 to 2.1%) |
| D, level 1 | 5,539 | 0.0535 | 0.1960 | 1544 of 5539 (26.7 to 29.1%) | 78 of 5539 (1.1 to 1.8%) |
<!-- /table:ar3c_map_four -->

<!-- table:ar3c_map_vs_D -->
| level | n | P_wrong | PR_wrong | D_wrong | fixed | broken | both | iv |
|---|---|---|---|---|---|---|---|---|
| level 0 | 4,877 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) | 17 of 4877 (0.22 to 0.56%) | 16 | 2 | 1 | do not overlap |
| level 0.1 | 4,488 | 10 of 4488 (0.12 to 0.41%) | 10 of 4488 (0.12 to 0.41%) | 23 of 4488 (0.34 to 0.77%) | 16 | 3 | 7 | overlap |
| level 0.2 | 5,041 | 23 of 5041 (0.30 to 0.68%) | 23 of 5041 (0.30 to 0.68%) | 38 of 5041 (0.5 to 1.0%) | 23 | 8 | 15 | overlap |
| level 0.5 | 4,625 | 39 of 4625 (0.6 to 1.2%) | 39 of 4625 (0.6 to 1.2%) | 47 of 4625 (0.8 to 1.3%) | 22 | 14 | 25 | overlap |
| level 1 | 5,539 | 95 of 5539 (1.4 to 2.1%) | 95 of 5539 (1.4 to 2.1%) | 78 of 5539 (1.1 to 1.8%) | 14 | 31 | 64 | overlap |
| all | 24,570 | 170 of 24570 (0.60 to 0.80%) | 170 of 24570 (0.60 to 0.80%) | 203 of 24570 (0.72 to 0.95%) | 91 | 58 | 112 | overlap |
<!-- /table:ar3c_map_vs_D -->

<!-- table:ar3c_reading0 -->
| level | n | ht | hc | rm | r5 | below | ct | cc |
|---|---|---|---|---|---|---|---|---|
| level 0 | 4,877 | 1.000 | 1.000 | 1.000 | 1.000 | 0 of 4877 (0.00 to 0.08%) | 0.073 | 0.073 |
| level 0.1 | 4,488 | 1.000 | 0.951 | 0.957 | 0.815 | 979 of 4488 (20.6 to 23.0%) | 0.074 | 0.118 |
| level 0.2 | 5,041 | 1.000 | 0.897 | 0.901 | 0.667 | 2498 of 5041 (48.2 to 50.9%) | 0.063 | 0.159 |
| level 0.5 | 4,625 | 1.000 | 0.763 | 0.767 | 0.451 | 3618 of 4625 (77.0 to 79.4%) | 0.068 | 0.285 |
| level 1 | 5,539 | 1.000 | 0.594 | 0.596 | 0.281 | 5168 of 5539 (92.6 to 93.9%) | 0.068 | 0.449 |
<!-- /table:ar3c_reading0 -->

<!-- table:ar3c_refine_all -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| start, painted voxels | 19,013 | 0.0406 | 0.1329 | 2240 of 19013 (11.3 to 12.2%) | 87 of 19013 (0.37 to 0.56%) |
| after R | 19,013 | 0.0406 | 0.1325 | 2230 of 19013 (11.3 to 12.2%) | 87 of 19013 (0.37 to 0.56%) |
| D (BFS) | 19,013 | 0.0440 | 0.1397 | 2601 of 19013 (13.2 to 14.2%) | 105 of 19013 (0.46 to 0.67%) |
<!-- /table:ar3c_refine_all -->

<!-- table:ar3c_refine_bins -->
| err | n | moved | closer | hurt |
|---|---|---|---|---|
| 0 to 0.05 | 11,423 | 3 of 11423 (0.01 to 0.08%) | 0 of 11423 (0.00 to 0.03%) | 0 of 11423 (0.00 to 0.03%) |
| 0.05 to 0.1 | 5,350 | 11 of 5350 (0.11 to 0.37%) | 1 of 5350 (0.00 to 0.11%) | 1 of 5350 (0.00 to 0.11%) |
| 0.1 to 0.2 | 2,042 | 20 of 2042 (0.6 to 1.5%) | 6 of 2042 (0.13 to 0.64%) | 1 of 2042 (0.01 to 0.28%) |
| 0.2 to 0.5 | 111 | 2 of 111 (0.5 to 6.3%) | 2 of 111 (0.5 to 6.3%) | 0 of 111 (0.0 to 3.3%) |
| 0.5 and above | 87 | 0 of 87 (0.0 to 4.2%) | 0 of 87 (0.0 to 4.2%) | 0 of 87 (0.0 to 4.2%) |
<!-- /table:ar3c_refine_bins -->

<!-- table:ar3c_refine_level -->
| level | n | closer | hurt | med_s | med_r | med_d | p95_s | p95_r | p95_d | w_s | w_r | w_d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| level 0 | 4,774 | 0 of 4774 (0.00 to 0.08%) | 0 of 4774 (0.00 to 0.08%) | 0.0161 | 0.0161 | 0.0170 | 0.0518 | 0.0518 | 0.0554 | 3 | 3 | 16 |
| level 0.1 | 4,152 | 0 of 4152 (0.00 to 0.09%) | 0 of 4152 (0.00 to 0.09%) | 0.0324 | 0.0324 | 0.0448 | 0.0704 | 0.0704 | 0.0968 | 9 | 9 | 20 |
| level 0.2 | 4,332 | 0 of 4332 (0.00 to 0.09%) | 1 of 4332 (0.00 to 0.13%) | 0.0433 | 0.0433 | 0.0559 | 0.0917 | 0.0917 | 0.1185 | 20 | 20 | 27 |
| level 0.5 | 3,196 | 2 of 3196 (0.02 to 0.23%) | 0 of 3196 (0.00 to 0.12%) | 0.0719 | 0.0718 | 0.0765 | 0.1440 | 0.1432 | 0.1647 | 25 | 25 | 26 |
| level 1 | 2,559 | 7 of 2559 (0.13 to 0.56%) | 1 of 2559 (0.01 to 0.22%) | 0.1022 | 0.1018 | 0.0723 | 0.1989 | 0.1986 | 0.1914 | 30 | 30 | 16 |
| all | 19,013 | 9 of 19013 (0.02 to 0.09%) | 2 of 19013 (0.00 to 0.04%) | 0.0406 | 0.0406 | 0.0440 | 0.1329 | 0.1325 | 0.1397 | 87 | 87 | 105 |
<!-- /table:ar3c_refine_level -->

<!-- table:ar3c_second -->
| level | n | moved | helped | hurt |
|---|---|---|---|---|
| level 0 | 1206 | 0 | 0 | 0 |
| level 0.1 | 1052 | 0 | 0 | 0 |
| level 0.2 | 1088 | 1 | 1 | 0 |
| level 0.5 | 754 | 1 | 0 | 0 |
| level 1 | 668 | 7 | 3 | 0 |
| all | 4768 | 9 | 4 | 0 |
<!-- /table:ar3c_second -->

<!-- table:ar3c_second_rule -->
| threshold | A_pred | A_helped | B_pred | B_right | B_base |
|---|---|---|---|---|---|
| 0.90 | 2471 | 2 | 2290 | 2288 of 2290 (99.7 to 100.0%) | 2 of 2292 (0.02 to 0.32%) |
<!-- /table:ar3c_second_rule -->

<!-- table:ar3c_second_wil -->
| what | v |
|---|---|
| moved by more than 0.05 deg | 9 of 4768 (0.10 to 0.36%) |
| brought closer by more than 0.05 deg | 4 of 4768 (0.03 to 0.22%) |
| pushed away by more than 0.05 deg | 0 of 4768 (0.00 to 0.08%) |
<!-- /table:ar3c_second_wil -->

<!-- table:ar3c_seed_cost -->
| what | n | per |
|---|---|---|
| entries from the seed path | 5,197 | 2,714 |
| entries from the grid search (fallback) | 172 | 339,891 |
<!-- /table:ar3c_seed_cost -->

<!-- table:ar3c_skip -->
| rule | threshold | A_share | B_skipped | B_evals | B_hurt | B_med | B_p95 | B_gt01 | B_wrong |
|---|---|---|---|---|---|---|---|---|---|
| start hit ratio >= t x entry best | 0.9 | 0.998 | 9299 of 9305 | 6,307,982 of 6,311,954 | 3 of 9305 (0.01 to 0.09%) | 0.0365 | 0.1378 | 1467 of 12184 (11.5 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
| start cost <= c | 0.3054 | 1.000 | 9305 of 9305 | 6,311,954 of 6,311,954 | 3 of 9305 (0.01 to 0.09%) | 0.0365 | 0.1378 | 1467 of 12184 (11.5 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
<!-- /table:ar3c_skip -->

<!-- table:ar3c_skip_level -->
| rule | level | n | skipped | hurt | skipped_evals |
|---|---|---|---|---|---|
| hit | level 0 | 2170 | 2170 | 0 | 1,451,098 |
| hit | level 0.1 | 1921 | 1920 | 0 | 1,287,104 |
| hit | level 0.2 | 2097 | 2097 | 0 | 1,410,390 |
| hit | level 0.5 | 1459 | 1456 | 1 | 990,108 |
| hit | level 1 | 1658 | 1656 | 2 | 1,169,282 |
| cost | level 0 | 2170 | 2170 | 0 | 1,451,098 |
| cost | level 0.1 | 1921 | 1921 | 0 | 1,287,766 |
| cost | level 0.2 | 2097 | 2097 | 0 | 1,410,390 |
| cost | level 0.5 | 1459 | 1459 | 1 | 992,094 |
| cost | level 1 | 1658 | 1658 | 2 | 1,170,606 |
<!-- /table:ar3c_skip_level -->

<!-- table:ar3c_sources -->
| group | n | grain | fit | seed | failed | none |
|---|---|---|---|---|---|---|
| all voxels | 24,570 | 19013 | 184 | 5369 | 4 | 0 |
| level 0 | 4,877 | 4774 | 0 | 103 | 0 | 0 |
| level 0.1 | 4,488 | 4152 | 2 | 334 | 0 | 0 |
| level 0.2 | 5,041 | 4332 | 1 | 707 | 1 | 0 |
| level 0.5 | 4,625 | 3196 | 48 | 1381 | 0 | 0 |
| level 1 | 5,539 | 2559 | 133 | 2844 | 3 | 0 |
<!-- /table:ar3c_sources -->

<!-- table:ar3c_start -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| all voxels | 24,570 | 0.0342 | 0.1304 | 2541 of 24570 (10.0 to 10.7%) | 170 of 24570 (0.60 to 0.80%) |
| level 0 | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| level 0.1 | 4,488 | 0.0313 | 0.0699 | 41 of 4488 (0.7 to 1.2%) | 10 of 4488 (0.12 to 0.41%) |
| level 0.2 | 5,041 | 0.0401 | 0.0896 | 151 of 5041 (2.6 to 3.5%) | 23 of 5041 (0.30 to 0.68%) |
| level 0.5 | 4,625 | 0.0571 | 0.1375 | 773 of 4625 (15.7 to 17.8%) | 39 of 4625 (0.6 to 1.2%) |
| level 1 | 5,539 | 0.0551 | 0.1875 | 1573 of 5539 (27.2 to 29.6%) | 95 of 5539 (1.4 to 2.1%) |
<!-- /table:ar3c_start -->

<!-- table:ar3c_stop -->
| rule | threshold | A_saved_share | B_saved | B_saved_share | B_hurt | B_med | B_p95 | B_gt01 | B_wrong |
|---|---|---|---|---|---|---|---|---|---|
| stop at the first failed block | none | 0.981 | 6,190,769 of 6,311,954 | 0.981 | 0 of 9305 (0.00 to 0.04%) | 0.0364 | 0.1376 | 1462 of 12184 (11.4 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
| stop when best cost <= c | 0.3054 | 0.981 | 6,190,989 of 6,311,954 | 0.981 | 0 of 9305 (0.00 to 0.04%) | 0.0364 | 0.1376 | 1462 of 12184 (11.4 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
<!-- /table:ar3c_stop -->

<!-- table:ar3c_stop_level -->
| rule | level | n | saved | total | hurt |
|---|---|---|---|---|---|
| first_fail | level 0 | 2170 | 1,422,888 | 1,451,098 | 0 |
| first_fail | level 0.1 | 1921 | 1,262,793 | 1,287,766 | 0 |
| first_fail | level 0.2 | 2097 | 1,383,085 | 1,410,390 | 0 |
| first_fail | level 0.5 | 1459 | 973,061 | 992,094 | 0 |
| first_fail | level 1 | 1658 | 1,148,942 | 1,170,606 | 0 |
| cost | level 0 | 2170 | 1,422,888 | 1,451,098 | 0 |
| cost | level 0.1 | 1921 | 1,262,793 | 1,287,766 | 0 |
| cost | level 0.2 | 2097 | 1,383,129 | 1,410,390 | 0 |
| cost | level 0.5 | 1459 | 973,127 | 992,094 | 0 |
| cost | level 1 | 1658 | 1,149,052 | 1,170,606 | 0 |
<!-- /table:ar3c_stop_level -->

<!-- table:ar3c_wall -->
| run | wall | label |
|---|---|---|
| P (paint, with set-up of pairing) | 7390 | contended: timing preflight refused (cpptools at 100% CPU, Python at 99% CPU); single run, one worker |
| R (one fit per painted voxel) | 841 | same process as P |
| R2 (second fit, every fourth voxel) | 219 | same process as P |
| D (BFS) | 9369 | contended: timing preflight refused (cpptools at 114% CPU, Python at 100% CPU); single run, one worker |
<!-- /table:ar3c_wall -->
