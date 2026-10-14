<!-- table:step7_candidates -->
| arm | images | median candidates | median passing the screen | median closest (deg) | closest within 1 deg | closest within 3 deg | median within 1 deg | median within 3 deg | median rank of closest | closest ranked first | closest in top quarter |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lenient rule | clean | 182 | 190 | 1.295 | 54 of 200 (21 to 34%) | 196 of 200 (95 to 99%) | 0 | 3 | 10 | 39 | 160 |
| lenient rule | no drop | 289 | 306 | 1.295 | 55 of 200 (22 to 34%) | 199 of 200 (97 to 100%) | 0 | 4 | 7 | 48 | 174 |
| C++ rule | clean | 125 | 128 | 1.310 | 52 of 200 (20 to 32%) | 194 of 200 (94 to 99%) | 0 | 3 | 13 | 27 | 147 |
| C++ rule | no drop | 194 | 204 | 1.310 | 54 of 200 (21 to 34%) | 196 of 200 (95 to 99%) | 0 | 3 | 14 | 40 | 161 |
| ray rule, 2-pixel window | clean | 181 | 189 | 1.295 | 54 of 200 (21 to 34%) | 196 of 200 (95 to 99%) | 0 | 3 | 10 | 39 | 160 |
| ray rule, 2-pixel window | no drop | 284 | 301 | 1.295 | 55 of 200 (22 to 34%) | 199 of 200 (97 to 100%) | 0 | 4 | 7 | 47 | 174 |
<!-- /table:step7_candidates -->

<!-- table:step7_candidates_paired -->
| arm | candidates: more / fewer / equal (no drop) | median ratio of candidates | closest: nearer / farther / same (no drop) | within 1 deg: both / only no drop / only clean / neither | within 3 deg: both / only no drop / only clean / neither | in top quarter: both / only no drop / only clean / neither |
|---|---|---|---|---|---|---|
| lenient rule | 200 / 0 / 0 | 1.64 | 15 / 11 / 174 | 54 / 1 / 0 / 145 | 196 / 3 / 0 / 1 | 155 / 19 / 5 / 21 |
| C++ rule | 199 / 0 / 1 | 1.59 | 12 / 7 / 181 | 52 / 2 / 0 / 146 | 194 / 2 / 0 / 4 | 140 / 21 / 7 / 32 |
| ray rule, 2-pixel window | 200 / 0 / 0 | 1.60 | 15 / 11 / 174 | 54 / 1 / 0 / 145 | 196 / 3 / 0 / 1 | 154 / 20 / 6 / 20 |
<!-- /table:step7_candidates_paired -->

<!-- table:step7_groups_off -->
| group | images | voxels | median candidates | median closest (deg) | closest within 1 deg | closest within 3 deg | closest in top quarter | median rank of closest | median rank / candidates | level-0 width (deg) | screen pass at 1 deg | screen pass at 3 deg |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fixed | clean | 31 | 176 | 1.784 | 1 | 29 | 22 | 23 | 0.127 | 0.235 | 0.980 | 0.358 |
| fixed | no drop | 31 | 246 | 1.621 | 1 | 31 | 23 | 16 | 0.056 | 0.321 | 0.984 | 0.459 |
| broken | clean | 12 | 180 | 1.610 | 0 | 12 | 6 | 47 | 0.260 | 0.236 | 0.982 | 0.363 |
| broken | no drop | 12 | 282 | 1.612 | 0 | 12 | 9 | 48 | 0.180 | 0.293 | 0.983 | 0.402 |
| both right | clean | 147 | 184 | 1.197 | 53 | 147 | 127 | 6 | 0.033 | 0.234 | 0.977 | 0.370 |
| both right | no drop | 147 | 306 | 1.196 | 54 | 147 | 138 | 4 | 0.012 | 0.314 | 0.985 | 0.437 |
| both wrong | clean | 10 | 162 | 2.344 | 0 | 8 | 5 | 44 | 0.298 | 0.226 | 0.980 | 0.340 |
| both wrong | no drop | 10 | 275 | 2.206 | 0 | 9 | 4 | 70 | 0.352 | 0.276 | 0.986 | 0.394 |
<!-- /table:step7_groups_off -->

<!-- table:step7_groups_on -->
| group | images | voxels | median candidates | median closest (deg) | closest within 1 deg | closest within 3 deg | closest in top quarter | median rank of closest | median rank / candidates | level-0 width (deg) | screen pass at 1 deg | screen pass at 3 deg |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fixed | clean | 32 | 126 | 1.513 | 1 | 31 | 19 | 23 | 0.210 | 0.171 | 0.978 | 0.264 |
| fixed | no drop | 32 | 208 | 1.481 | 1 | 32 | 25 | 18 | 0.086 | 0.210 | 0.988 | 0.338 |
| broken | clean | 14 | 124 | 1.493 | 1 | 14 | 8 | 28 | 0.219 | 0.174 | 0.969 | 0.240 |
| broken | no drop | 14 | 222 | 1.664 | 1 | 14 | 8 | 46 | 0.202 | 0.212 | 0.973 | 0.270 |
| both right | clean | 130 | 128 | 1.175 | 50 | 129 | 107 | 8 | 0.063 | 0.171 | 0.966 | 0.258 |
| both right | no drop | 130 | 192 | 1.168 | 52 | 130 | 114 | 7 | 0.043 | 0.209 | 0.978 | 0.299 |
| both wrong | clean | 24 | 116 | 2.067 | 0 | 20 | 13 | 26 | 0.209 | 0.168 | 0.968 | 0.245 |
| both wrong | no drop | 24 | 172 | 2.067 | 0 | 20 | 14 | 40 | 0.225 | 0.202 | 0.980 | 0.285 |
<!-- /table:step7_groups_on -->

<!-- table:step7_groups_paired -->
| arm | group | voxels | in top quarter: both / only no drop / only clean / neither | within 3 deg: both / only no drop / only clean / neither | level-0 width, clean (deg) | level-0 width, no drop (deg) | width ratio | ratio 95% interval (voxel bootstrap) | screen pass at 3 deg, clean | screen pass at 3 deg, no drop |
|---|---|---|---|---|---|---|---|---|---|---|
| lenient rule | fixed | 31 | 20 / 3 / 2 / 6 | 29 / 2 / 0 / 0 | 0.235 | 0.321 | 1.36 | 1.31 to 1.43 | 0.358 | 0.459 |
| lenient rule | broken | 12 | 6 / 3 / 0 / 3 | 12 / 0 / 0 / 0 | 0.236 | 0.293 | 1.24 | 1.17 to 1.35 | 0.363 | 0.402 |
| lenient rule | both right | 147 | 125 / 13 / 2 / 7 | 147 / 0 / 0 / 0 | 0.234 | 0.314 | 1.34 | 1.31 to 1.37 | 0.370 | 0.437 |
| lenient rule | both wrong | 10 | 4 / 0 / 1 / 5 | 8 / 1 / 0 / 1 | 0.226 | 0.276 | 1.22 | 1.18 to 1.32 | 0.340 | 0.394 |
| C++ rule | fixed | 32 | 19 / 6 / 0 / 7 | 31 / 1 / 0 / 0 | 0.171 | 0.210 | 1.23 | 1.19 to 1.29 | 0.264 | 0.338 |
| C++ rule | broken | 14 | 6 / 2 / 2 / 4 | 14 / 0 / 0 / 0 | 0.174 | 0.212 | 1.22 | 1.14 to 1.34 | 0.240 | 0.270 |
| C++ rule | both right | 130 | 103 / 11 / 4 / 12 | 129 / 1 / 0 / 0 | 0.171 | 0.209 | 1.22 | 1.20 to 1.25 | 0.258 | 0.299 |
| C++ rule | both wrong | 24 | 12 / 2 / 1 / 9 | 20 / 0 / 0 / 4 | 0.168 | 0.202 | 1.20 | 1.16 to 1.24 | 0.245 | 0.285 |
| ray rule, 2-pixel window | fixed | 35 | 18 / 5 / 4 / 8 | 32 / 3 / 0 / 0 | 0.234 | 0.317 | 1.35 | 1.30 to 1.41 | 0.336 | 0.411 |
| ray rule, 2-pixel window | broken | 9 | 5 / 2 / 0 / 2 | 9 / 0 / 0 / 0 | 0.232 | 0.300 | 1.29 | 1.18 to 1.40 | 0.356 | 0.396 |
| ray rule, 2-pixel window | both right | 151 | 129 / 12 / 2 / 8 | 151 / 0 / 0 / 0 | 0.234 | 0.310 | 1.32 | 1.30 to 1.35 | 0.373 | 0.426 |
| ray rule, 2-pixel window | both wrong | 5 | 2 / 1 / 0 / 2 | 4 / 0 / 0 / 1 | 0.235 | 0.287 | 1.22 | 1.14 to 1.36 | 0.320 | 0.360 |
<!-- /table:step7_groups_paired -->

<!-- table:step7_groups_rayR2 -->
| group | images | voxels | median candidates | median closest (deg) | closest within 1 deg | closest within 3 deg | closest in top quarter | median rank of closest | median rank / candidates | level-0 width (deg) | screen pass at 1 deg | screen pass at 3 deg |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fixed | clean | 35 | 168 | 1.833 | 1 | 32 | 22 | 39 | 0.209 | 0.234 | 0.978 | 0.336 |
| fixed | no drop | 35 | 262 | 1.713 | 1 | 35 | 23 | 39 | 0.164 | 0.317 | 0.986 | 0.411 |
| broken | clean | 9 | 183 | 1.665 | 0 | 9 | 5 | 43 | 0.199 | 0.232 | 0.960 | 0.356 |
| broken | no drop | 9 | 280 | 1.681 | 0 | 9 | 7 | 36 | 0.156 | 0.300 | 0.969 | 0.396 |
| both right | clean | 151 | 184 | 1.197 | 53 | 151 | 131 | 7 | 0.039 | 0.234 | 0.979 | 0.373 |
| both right | no drop | 151 | 294 | 1.196 | 54 | 151 | 141 | 4 | 0.014 | 0.310 | 0.986 | 0.426 |
| both wrong | clean | 5 | 159 | 1.916 | 0 | 4 | 2 | 76 | 0.478 | 0.235 | 0.980 | 0.320 |
| both wrong | no drop | 5 | 261 | 1.916 | 0 | 4 | 3 | 73 | 0.222 | 0.287 | 0.980 | 0.360 |
<!-- /table:step7_groups_rayR2 -->

<!-- table:step7_rise_ratio -->
| arm | cost | 0.03 deg | 0.1 deg | 0.3 deg | 0.5 deg | 1 deg | 2 deg | 3 deg |
|---|---|---|---|---|---|---|---|---|
| lenient rule | level-0 cost | 0.930 | 0.746 | 0.806 | 0.895 | 0.976 | 0.999 | 1.000 |
| lenient rule | final cost | 0.836 | 0.838 | 0.899 | 0.925 | 0.944 | 0.949 | 0.951 |
| C++ rule | level-0 cost | 0.708 | 0.730 | 0.873 | 0.937 | 0.985 | 1.000 | 1.000 |
| C++ rule | final cost | 0.830 | 0.830 | 0.899 | 0.926 | 0.944 | 0.951 | 0.951 |
| ray rule, 2-pixel window | level-0 cost | 0.930 | 0.754 | 0.808 | 0.899 | 0.979 | 0.999 | 1.000 |
| ray rule, 2-pixel window | final cost | 0.830 | 0.830 | 0.899 | 0.926 | 0.944 | 0.951 | 0.951 |
<!-- /table:step7_rise_ratio -->

<!-- table:step7_screen_pass -->
| arm | images | truth | 0.03 deg | 0.1 deg | 0.3 deg | 0.5 deg | 1 deg | 2 deg | 3 deg |
|---|---|---|---|---|---|---|---|---|---|
| lenient rule | clean | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.978 | 0.692 | 0.366 |
| lenient rule | no drop | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.985 | 0.752 | 0.436 |
| C++ rule | clean | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.968 | 0.545 | 0.256 |
| C++ rule | no drop | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.979 | 0.617 | 0.301 |
| ray rule, 2-pixel window | clean | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.978 | 0.690 | 0.364 |
| ray rule, 2-pixel window | no drop | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.985 | 0.746 | 0.421 |
<!-- /table:step7_screen_pass -->

<!-- table:step7_width -->
| arm | cost | images | truth | 0.03 deg | 0.1 deg | 0.3 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | width, grid (deg) | width, interpolated (deg) | median per-voxel width, grid (deg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lenient rule | level-0 cost | clean | 0.0000 | 0.0303 | 0.1538 | 0.6000 | 0.8154 | 0.9545 | 0.9859 | 1.0000 | 0.3 | 0.235 | 0.3 |
| lenient rule | level-0 cost | no drop | 0.0000 | 0.0282 | 0.1148 | 0.4839 | 0.7302 | 0.9318 | 0.9848 | 1.0000 | 0.5 | 0.310 | 0.5 |
| lenient rule | final cost | clean | 0.0662 | 0.1920 | 0.4887 | 0.8432 | 0.9318 | 0.9833 | 0.9995 | 1.0000 | 0.3 | 0.115 | 0.3 |
| lenient rule | final cost | no drop | 0.1120 | 0.2172 | 0.4660 | 0.8103 | 0.9128 | 0.9778 | 0.9978 | 1.0000 | 0.3 | 0.133 | 0.3 |
| C++ rule | level-0 cost | clean | 0.0000 | 0.0435 | 0.2656 | 0.7463 | 0.8923 | 0.9701 | 1.0000 | 1.0000 | 0.3 | 0.171 | 0.3 |
| C++ rule | level-0 cost | no drop | 0.0000 | 0.0308 | 0.1940 | 0.6515 | 0.8361 | 0.9552 | 1.0000 | 1.0000 | 0.3 | 0.208 | 0.3 |
| C++ rule | final cost | clean | 0.0662 | 0.1947 | 0.5081 | 0.8589 | 0.9392 | 0.9852 | 1.0000 | 1.0000 | 0.3 | 0.108 | 0.3 |
| C++ rule | final cost | no drop | 0.1120 | 0.2186 | 0.4790 | 0.8247 | 0.9205 | 0.9799 | 1.0000 | 1.0000 | 0.3 | 0.128 | 0.3 |
| ray rule, 2-pixel window | level-0 cost | clean | 0.0000 | 0.0303 | 0.1538 | 0.6000 | 0.8154 | 0.9545 | 0.9859 | 1.0000 | 0.3 | 0.235 | 0.3 |
| ray rule, 2-pixel window | level-0 cost | no drop | 0.0000 | 0.0282 | 0.1159 | 0.4848 | 0.7333 | 0.9344 | 0.9848 | 1.0000 | 0.5 | 0.309 | 0.5 |
| ray rule, 2-pixel window | final cost | clean | 0.0662 | 0.1947 | 0.5081 | 0.8589 | 0.9392 | 0.9852 | 1.0000 | 1.0000 | 0.3 | 0.108 | 0.3 |
| ray rule, 2-pixel window | final cost | no drop | 0.1120 | 0.2186 | 0.4790 | 0.8247 | 0.9205 | 0.9799 | 1.0000 | 1.0000 | 0.3 | 0.128 | 0.3 |
<!-- /table:step7_width -->

<!-- table:step7_width_mean -->
| arm | images | truth | 0.03 deg | 0.1 deg | 0.3 deg | 0.5 deg | 1 deg | 2 deg | 3 deg |
|---|---|---|---|---|---|---|---|---|---|
| lenient rule | clean | 0.0051 | 0.0305 | 0.1567 | 0.5804 | 0.7967 | 0.9523 | 0.9890 | 0.9951 |
| lenient rule | no drop | 0.0051 | 0.0259 | 0.1159 | 0.4655 | 0.7089 | 0.9276 | 0.9835 | 0.9926 |
| C++ rule | clean | 0.0051 | 0.0422 | 0.2578 | 0.7170 | 0.8731 | 0.9700 | 0.9929 | 0.9969 |
| C++ rule | no drop | 0.0051 | 0.0344 | 0.1928 | 0.6197 | 0.8129 | 0.9544 | 0.9895 | 0.9954 |
| ray rule, 2-pixel window | clean | 0.0051 | 0.0306 | 0.1569 | 0.5806 | 0.7968 | 0.9524 | 0.9890 | 0.9951 |
| ray rule, 2-pixel window | no drop | 0.0051 | 0.0260 | 0.1162 | 0.4682 | 0.7115 | 0.9287 | 0.9839 | 0.9929 |
<!-- /table:step7_width_mean -->

<!-- table:step7_width_paired -->
| arm | cost | voxels | no drop wider | no drop narrower | same | median per-voxel width, clean (deg) | median per-voxel width, no drop (deg) | pooled width, clean (deg) | pooled width, no drop (deg) | pooled width ratio (no drop / clean) | ratio 95% interval (voxel bootstrap) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lenient rule | level-0 cost | 200 | 200 | 0 | 0 | 0.235 | 0.306 | 0.235 | 0.310 | 1.323 | 1.31 to 1.36 |
| lenient rule | final cost | 200 | 200 | 0 | 0 | 0.115 | 0.130 | 0.115 | 0.133 | 1.161 | 1.15 to 1.18 |
| C++ rule | level-0 cost | 200 | 200 | 0 | 0 | 0.172 | 0.204 | 0.171 | 0.208 | 1.220 | 1.20 to 1.24 |
| C++ rule | final cost | 200 | 200 | 0 | 0 | 0.108 | 0.124 | 0.108 | 0.128 | 1.181 | 1.17 to 1.20 |
| ray rule, 2-pixel window | level-0 cost | 200 | 200 | 0 | 0 | 0.235 | 0.305 | 0.235 | 0.309 | 1.320 | 1.30 to 1.35 |
| ray rule, 2-pixel window | final cost | 200 | 200 | 0 | 0 | 0.108 | 0.124 | 0.108 | 0.128 | 1.181 | 1.17 to 1.20 |
<!-- /table:step7_width_paired -->
