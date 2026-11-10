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
