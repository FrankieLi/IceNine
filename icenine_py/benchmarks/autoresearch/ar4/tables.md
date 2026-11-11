<!-- table:ar4_boundary -->
| score | auc | wrong | flagged |
|---|---|---|---|
| B0 (hit ratio) | 0.974 | 31 of 31 | 1735 of 17859 |
| B1 (logistic, summary numbers) | 0.882 | 31 of 31 | 9487 of 17859 |
| N1 (set network) | 0.991 | 23 of 31 | 189 of 17859 |
| share of peaks with a lit pixel | 0.496 | 31 of 31 | 17859 of 17859 |
<!-- /table:ar4_boundary -->

<!-- table:ar4_cuts -->
| score | cut | flagged | wrong | share |
|---|---|---|---|---|
| B0 (hit ratio) | 0.926 | 1741 | 83 | 3.54% |
| B1 (logistic, summary numbers) | 0.976 | 25009 | 83 | 50.89% |
| N1 (set network) | 0.749 | 110 | 83 | 0.22% |
| share of peaks with a lit pixel | 1.000 | 49140 | 87 | 100.00% |
<!-- /table:ar4_cuts -->

<!-- table:ar4_lit -->
| sample | all_lit | wrong_all_lit | lowest | median |
|---|---|---|---|---|
| big | 77.1% | 19 of 24 | 0.972 | 1.000 |
| grad | 69.8% | 48 of 63 | 0.900 | 1.000 |
| mid | 79.1% | 24 of 31 | 0.984 | 1.000 |
<!-- /table:ar4_lit -->

<!-- table:ar4_missed -->
| voxel | err | hit | n1 | bnd |
|---|---|---|---|---|
| 97 | 55.06 | 0.903 | 0.986 | yes |
| 2995 | 34.40 | 0.908 | 0.959 | yes |
| 3302 | 54.20 | 0.889 | 0.912 | yes |
| 4758 | 52.10 | 0.904 | 0.982 | yes |
| 5807 | 38.94 | 0.905 | 0.828 | yes |
| 6419 | 52.50 | 0.910 | 0.899 | yes |
| 8427 | 55.99 | 0.899 | 0.973 | yes |
| 21630 | 54.73 | 0.900 | 0.926 | yes |
<!-- /table:ar4_missed -->

<!-- table:ar4_pairs -->
| sample | pairs | after | before | try_start | seed_end | seed_cand | turn_01 | turn_03 | turn_05 |
|---|---|---|---|---|---|---|---|---|---|
| big | 49103 | 24570 | 1056 | 18665 | 465 | 660 | 1229 | 1229 | 1229 |
| grad | 94306 | 24570 | 752 | 60381 | 2356 | 2560 | 1229 | 1229 | 1229 |
| mid | 86242 | 24570 | 1429 | 50099 | 3217 | 3240 | 1229 | 1229 | 1229 |
<!-- /table:ar4_pairs -->

<!-- table:ar4_use1_big -->
| score | wrong | flagged | share | auc | rule |
|---|---|---|---|---|---|
| B0 (hit ratio) | 21 of 24 (0.69 to 0.96) | 553 of 24570 | 2.25% | 0.985 | no |
| B1 (logistic, summary numbers) | 22 of 24 (0.74 to 0.98) | 12800 of 24570 | 52.10% | 0.842 | no |
| N1 (set network) | 22 of 24 (0.74 to 0.98) | 25 of 24570 | 0.10% | 1.000 | yes |
| share of peaks with a lit pixel | 24 of 24 (0.86 to 1.00) | 24570 of 24570 | 100.00% | 0.482 | no |
<!-- /table:ar4_use1_big -->

<!-- table:ar4_use1_grad -->
| score | wrong | flagged | share | auc | rule |
|---|---|---|---|---|---|
| B0 (hit ratio) | 62 of 63 (0.92 to 1.00) | 1188 of 24570 | 4.84% | 0.987 | no |
| B1 (logistic, summary numbers) | 61 of 63 (0.89 to 0.99) | 12209 of 24570 | 49.69% | 0.918 | no |
| N1 (set network) | 61 of 63 (0.89 to 0.99) | 85 of 24570 | 0.35% | 1.000 | no |
| share of peaks with a lit pixel | 63 of 63 (0.94 to 1.00) | 24570 of 24570 | 100.00% | 0.478 | no |
<!-- /table:ar4_use1_grad -->

<!-- table:ar4_use1_mid -->
| score | wrong | flagged | share | auc | rule |
|---|---|---|---|---|---|
| B0 (hit ratio) | 31 of 31 (0.89 to 1.00) | 1760 of 24570 | 7.16% | 0.981 | no |
| B1 (logistic, summary numbers) | 31 of 31 (0.89 to 1.00) | 13653 of 24570 | 55.57% | 0.875 | no |
| N1 (set network) | 23 of 31 (0.57 to 0.86) | 190 of 24570 | 0.77% | 0.994 | no |
| share of peaks with a lit pixel | 31 of 31 (0.89 to 1.00) | 24570 of 24570 | 100.00% | 0.506 | no |
<!-- /table:ar4_use1_mid -->

<!-- table:ar4_use2_grad -->
| score | right | wrong | saved | refused |
|---|---|---|---|---|
| B0 (hit ratio) | 0 of 13 (0.00 to 0.23) | 0 of 3 (0.00 to 0.56) | 0 | 116 of 1263 (94 right) |
| B1 (logistic, summary numbers) | 6 of 13 (0.23 to 0.71) | 0 of 3 (0.00 to 0.56) | 2,061,839 | 356 of 1263 (337 right) |
| N1 (set network) | 11 of 13 (0.58 to 0.96) | 0 of 3 (0.00 to 0.56) | 3,787,006 | 31 of 1263 (10 right) |
| share of peaks with a lit pixel | 4 of 13 (0.13 to 0.58) | 1 of 3 (0.06 to 0.79) | 1,380,986 | 347 of 1263 (340 right) |
<!-- /table:ar4_use2_grad -->

<!-- table:ar4_use2_mid -->
| score | right | wrong | saved | refused |
|---|---|---|---|---|
| B0 (hit ratio) | 0 of 10 (0.00 to 0.28) | 0 of 5 (0.00 to 0.43) | 0 | 141 of 1596 (130 right) |
| B1 (logistic, summary numbers) | 3 of 10 (0.11 to 0.60) | 0 of 5 (0.00 to 0.43) | 1,017,386 | 582 of 1596 (571 right) |
| N1 (set network) | 7 of 10 (0.40 to 0.89) | 0 of 5 (0.00 to 0.43) | 2,365,778 | 29 of 1596 (20 right) |
| share of peaks with a lit pixel | 2 of 10 (0.06 to 0.51) | 4 of 5 (0.38 to 0.96) | 665,115 | 476 of 1596 (472 right) |
<!-- /table:ar4_use2_mid -->

<!-- table:ar4_whatif -->
| score | p95 | p100 |
|---|---|---|
| B0 (hit ratio) | 1016 (4.14%) | 1029 (4.19%) |
| B1 (logistic, summary numbers) | 9548 (38.86%) | 11910 (48.47%) |
| N1 (set network) | 779 (3.17%) | 854 (3.48%) |
| share of peaks with a lit pixel | 24570 (100.00%) | 24570 (100.00%) |
<!-- /table:ar4_whatif -->
