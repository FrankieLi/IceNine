<!-- table:q8_sanity_cost -->
| variant | 6 voxels, lenient (sanity_checks) | lenient | C++ rule | ray 2px | ray geometric |
|---|---|---|---|---|---|
| clean | 0.0674 | 0.0643 | 0.0643 | 0.0643 | 0.0643 |
| realistic | 0.2586 | 0.2646 | 0.2831 | 0.2831 | 0.2825 |
| realistic_nodrop | 0.1168 | 0.1162 | 0.1162 | 0.1162 | 0.1162 |
<!-- /table:q8_sanity_cost -->

<!-- table:q8_sanity_rejected -->
| variant | setting | det0_overlaps | C++ rule | C++ rule 95% interval | ray 2px | ray 2px 95% interval | ray geometric | ray geometric 95% interval |
|---|---|---|---|---|---|---|---|---|
| clean | screen | 524 | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% |
| clean | local | 521 | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% |
| realistic | screen | 475 | 43 of 475 (9.1%) | 6.79 to 11.97% | 43 of 475 (9.1%) | 6.79 to 11.97% | 39 of 475 (8.2%) | 6.06 to 11.03% |
| realistic | local | 471 | 41 of 471 (8.7%) | 6.48 to 11.60% | 40 of 471 (8.5%) | 6.30 to 11.36% | 36 of 471 (7.6%) | 5.57 to 10.40% |
| realistic_nodrop | screen | 524 | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% |
| realistic_nodrop | local | 521 | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% |
<!-- /table:q8_sanity_rejected -->
<!-- table:q8_e0_truth_cost -->
| msin | images | lenient | C++ rule | ray 2px | ray geometric |
|---|---|---|---|---|---|
| 0.3 | clean | 0.0682 | 0.0682 | 0.0682 | 0.0682 |
| 0.3 | realistic (10% drop) | 0.2627 | 0.2794 | 0.2793 | 0.2793 |
| 0.3 | realistic, no drop | 0.1148 | 0.1148 | 0.1148 | 0.1148 |
| 0.0 | clean | 0.0662 | 0.0662 | 0.0662 | 0.0662 |
| 0.0 | realistic (10% drop) | 0.2531 | 0.2681 | 0.2677 | 0.2677 |
| 0.0 | realistic, no drop | 0.1120 | 0.1120 | 0.1120 | 0.1120 |
<!-- /table:q8_e0_truth_cost -->

<!-- table:q8_e0_truth_rejected -->
| msin | images | det0_overlaps | C++ rule | C++ rule 95% interval | ray 2px | ray 2px 95% interval | ray geometric | ray geometric 95% interval |
|---|---|---|---|---|---|---|---|---|
| 0.3 | clean | 12439 | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% |
| 0.3 | realistic (10% drop) | 11174 | 974 of 11174 (8.7%) | 8.21 to 9.25% | 965 of 11174 (8.6%) | 8.13 to 9.17% | 964 of 11174 (8.6%) | 8.12 to 9.16% |
| 0.3 | realistic, no drop | 12439 | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% |
| 0 | clean | 12932 | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% |
| 0 | realistic (10% drop) | 11667 | 975 of 11667 (8.4%) | 7.87 to 8.87% | 966 of 11667 (8.3%) | 7.79 to 8.79% | 965 of 11667 (8.3%) | 7.79 to 8.78% |
| 0 | realistic, no drop | 12932 | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% |
<!-- /table:q8_e0_truth_rejected -->
<!-- table:q8_e0_paired_rules -->
| images | rule | reference | wrong_rule | wrong_ref | fixed | broken | both_wrong |
|---|---|---|---|---|---|---|---|
| clean | C++ rule, min_sin_eta 0.3 | lenient rule, min_sin_eta 0.3 | 61 | 35 | 15 | 41 | 20 |
| realistic (10% drop) | C++ rule, min_sin_eta 0.3 | lenient rule, min_sin_eta 0.3 | 42 | 27 | 10 | 25 | 17 |
| realistic, no drop | C++ rule, min_sin_eta 0.3 | lenient rule, min_sin_eta 0.3 | 39 | 16 | 5 | 28 | 11 |
| clean | C++ rule, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 56 | 41 | 15 | 30 | 26 |
| realistic (10% drop) | C++ rule, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 56 | 25 | 7 | 38 | 18 |
| realistic, no drop | C++ rule, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 38 | 22 | 10 | 26 | 12 |
| clean | ray rule 2-pixel window, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 40 | 41 | 10 | 9 | 31 |
| realistic (10% drop) | ray rule 2-pixel window, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 32 | 25 | 10 | 17 | 15 |
| realistic, no drop | ray rule 2-pixel window, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 14 | 22 | 15 | 7 | 7 |
| clean | ray rule geometric window, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 41 | 41 | 11 | 11 | 30 |
| realistic (10% drop) | ray rule geometric window, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 29 | 25 | 14 | 18 | 11 |
| realistic, no drop | ray rule geometric window, min_sin_eta 0 | lenient rule, min_sin_eta 0 | 14 | 22 | 14 | 6 | 8 |
| clean | ray rule 2-pixel window, min_sin_eta 0 | C++ rule, min_sin_eta 0 | 40 | 56 | 32 | 16 | 24 |
| realistic (10% drop) | ray rule 2-pixel window, min_sin_eta 0 | C++ rule, min_sin_eta 0 | 32 | 56 | 38 | 14 | 18 |
| realistic, no drop | ray rule 2-pixel window, min_sin_eta 0 | C++ rule, min_sin_eta 0 | 14 | 38 | 30 | 6 | 8 |
| clean | ray rule geometric window, min_sin_eta 0 | C++ rule, min_sin_eta 0 | 41 | 56 | 31 | 16 | 25 |
| realistic (10% drop) | ray rule geometric window, min_sin_eta 0 | C++ rule, min_sin_eta 0 | 29 | 56 | 37 | 10 | 19 |
| realistic, no drop | ray rule geometric window, min_sin_eta 0 | C++ rule, min_sin_eta 0 | 14 | 38 | 29 | 5 | 9 |
| clean | ray rule geometric window, min_sin_eta 0 | ray rule 2-pixel window, min_sin_eta 0 | 41 | 40 | 4 | 5 | 36 |
| realistic (10% drop) | ray rule geometric window, min_sin_eta 0 | ray rule 2-pixel window, min_sin_eta 0 | 29 | 32 | 12 | 9 | 20 |
| realistic, no drop | ray rule geometric window, min_sin_eta 0 | ray rule 2-pixel window, min_sin_eta 0 | 14 | 14 | 3 | 3 | 11 |
<!-- /table:q8_e0_paired_rules -->

<!-- table:q8_e0_paired_variants -->
| arm | comparison | fixed | broken | both_wrong | both_right |
|---|---|---|---|---|---|
| lenient rule, min_sin_eta 0.3 | no drop vs clean | 28 | 9 | 7 | 156 |
| lenient rule, min_sin_eta 0.3 | no drop vs old realistic | 17 | 6 | 10 | 167 |
| lenient rule, min_sin_eta 0.3 | old realistic vs clean | 24 | 16 | 11 | 149 |
| C++ rule, min_sin_eta 0.3 | no drop vs clean | 37 | 15 | 24 | 124 |
| C++ rule, min_sin_eta 0.3 | no drop vs old realistic | 21 | 18 | 21 | 140 |
| C++ rule, min_sin_eta 0.3 | old realistic vs clean | 34 | 15 | 27 | 124 |
| lenient rule, min_sin_eta 0 | no drop vs clean | 31 | 12 | 10 | 147 |
| lenient rule, min_sin_eta 0 | no drop vs old realistic | 14 | 11 | 11 | 164 |
| lenient rule, min_sin_eta 0 | old realistic vs clean | 27 | 11 | 14 | 148 |
| C++ rule, min_sin_eta 0 | no drop vs clean | 32 | 14 | 24 | 130 |
| C++ rule, min_sin_eta 0 | no drop vs old realistic | 33 | 15 | 23 | 129 |
| C++ rule, min_sin_eta 0 | old realistic vs clean | 27 | 27 | 29 | 117 |
| ray rule 2-pixel window, min_sin_eta 0 | no drop vs clean | 35 | 9 | 5 | 151 |
| ray rule 2-pixel window, min_sin_eta 0 | no drop vs old realistic | 26 | 8 | 6 | 160 |
| ray rule 2-pixel window, min_sin_eta 0 | old realistic vs clean | 23 | 15 | 17 | 145 |
| ray rule geometric window, min_sin_eta 0 | no drop vs clean | 35 | 8 | 6 | 151 |
| ray rule geometric window, min_sin_eta 0 | no drop vs old realistic | 24 | 9 | 5 | 162 |
| ray rule geometric window, min_sin_eta 0 | old realistic vs clean | 30 | 18 | 11 | 141 |
<!-- /table:q8_e0_paired_variants -->

<!-- table:q8_e0_wrong -->
| arm | images | wrong | interval | evals | runtime | wall_min | med_err | stored_earlier |
|---|---|---|---|---|---|---|---|---|
| lenient rule, min_sin_eta 0.3 | clean | 35 of 200 | 0.129 to 0.234 | 48623 | 24.9 | 8.5 | 0.0123 | 40 |
| lenient rule, min_sin_eta 0.3 | realistic (10% drop) | 27 of 200 | 0.094 to 0.189 | 51454 | 26.9 | 9.3 | 0.0081 | 25 |
| lenient rule, min_sin_eta 0.3 | realistic, no drop | 16 of 200 | 0.050 to 0.126 | 51825 | 27.1 | 9.3 | 0.0116 | n/a |
| C++ rule, min_sin_eta 0.3 | clean | 61 of 200 | 0.245 to 0.372 | 47174 | 23.6 | 8.1 | 0.0193 | 61 |
| C++ rule, min_sin_eta 0.3 | realistic (10% drop) | 42 of 200 | 0.159 to 0.272 | 48731 | 24.0 | 8.3 | 0.0148 | 42 |
| C++ rule, min_sin_eta 0.3 | realistic, no drop | 39 of 200 | 0.146 to 0.255 | 49227 | 24.3 | 8.4 | 0.0121 | n/a |
| lenient rule, min_sin_eta 0 | clean | 41 of 200 | 0.155 to 0.266 | 48700 | 23.5 | 8.1 | 0.0130 | 39 |
| lenient rule, min_sin_eta 0 | realistic (10% drop) | 25 of 200 | 0.086 to 0.178 | 51520 | 25.0 | 8.6 | 0.0100 | 18 |
| lenient rule, min_sin_eta 0 | realistic, no drop | 22 of 200 | 0.074 to 0.161 | 51901 | 25.5 | 8.8 | 0.0080 | n/a |
| C++ rule, min_sin_eta 0 | clean | 56 of 200 | 0.222 to 0.346 | 47230 | 23.5 | 8.1 | 0.0223 | 56 |
| C++ rule, min_sin_eta 0 | realistic (10% drop) | 56 of 200 | 0.222 to 0.346 | 48769 | 24.4 | 8.4 | 0.0132 | 56 |
| C++ rule, min_sin_eta 0 | realistic, no drop | 38 of 200 | 0.142 to 0.250 | 49281 | 24.4 | 8.4 | 0.0109 | n/a |
| ray rule 2-pixel window, min_sin_eta 0 | clean | 40 of 200 | 0.150 to 0.261 | 48641 | 26.5 | 9.1 | 0.0097 | 40 |
| ray rule 2-pixel window, min_sin_eta 0 | realistic (10% drop) | 32 of 200 | 0.116 to 0.217 | 50915 | 27.2 | 9.4 | 0.0087 | 32 |
| ray rule 2-pixel window, min_sin_eta 0 | realistic, no drop | 14 of 200 | 0.042 to 0.114 | 51652 | 27.0 | 9.3 | 0.0072 | n/a |
| ray rule geometric window, min_sin_eta 0 | clean | 41 of 200 | 0.155 to 0.266 | 48649 | 25.1 | 8.6 | 0.0114 | 41 |
| ray rule geometric window, min_sin_eta 0 | realistic (10% drop) | 29 of 200 | 0.103 to 0.200 | 50967 | 26.7 | 9.2 | 0.0100 | 29 |
| ray rule geometric window, min_sin_eta 0 | realistic, no drop | 14 of 200 | 0.042 to 0.114 | 51677 | 27.1 | 9.3 | 0.0086 | n/a |
<!-- /table:q8_e0_wrong -->
