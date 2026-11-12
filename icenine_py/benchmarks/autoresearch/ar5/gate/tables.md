<!-- table:curve -->
| h | mid S11n | grad S11n | big S11 | grad S11 | mid S11 | held S11v2 |
|---|---|---|---|---|---|---|
| 0.90 | 94.0% / 38 | 98.0% / 119 | 98.2% / 52 | 98.1% / 112 | 94.5% / 37 | 93.3% / 27 |
| 0.91 | 89.3% / 19 | 95.4% / 36 | 94.1% / 36 | 95.5% / 35 | 89.6% / 18 | 89.4% / 11 |
| 0.92 | 84.6% / 6 | 92.9% / 8 | 90.0% / 14 | 93.1% / 9 | 84.8% / 6 | 85.5% / 1 |
| 0.93 | 79.4% / 3 | 90.1% / 2 | 85.9% / 4 | 90.5% / 2 | 80.0% / 3 | 81.6% / 0 |
| 0.94 | 74.4% / 1 | 87.0% / 0 | 81.1% / 1 | 87.2% / 0 | 74.8% / 1 | 77.3% / 0 |
| 0.95 | 68.9% / 0 | 83.3% / 0 | 75.9% / 0 | 83.5% / 0 | 69.2% / 0 | 72.3% / 0 |
| 0.96 | 61.8% / 0 | 77.9% / 0 | 69.3% / 0 | 78.2% / 0 | 62.2% / 0 | 65.3% / 0 |
| 0.97 | 53.7% / 0 | 70.8% / 0 | 61.7% / 0 | 70.2% / 0 | 54.3% / 0 | 56.9% / 0 |
| 0.98 | 42.5% / 0 | 58.1% / 0 | 50.9% / 0 | 58.1% / 0 | 43.1% / 0 | 44.8% / 0 |
| 0.99 | 28.5% / 0 | 39.4% / 0 | 36.0% / 0 | 39.3% / 0 | 28.8% / 0 | 29.5% / 0 |
<!-- /table:curve -->

<!-- table:fixed_hits -->
| sample | run | voxels the full pass fixed | lowest own hit ratio | highest own hit ratio | of them with own hit ratio at least 0.95 | highest own hit ratio of a voxel wrong before / after the pass |
|---|---|---|---|---|---|---|
| mid | S11n | 43 | 0.889 | 0.944 | 0 | 0.944 / 0.910 |
| first gradient | S11n | 159 | 0.841 | 0.936 | 0 | 0.936 / 0.925 |
| big | S11 | 53 | 0.899 | 0.944 | 0 | 0.944 / 0.933 |
| first gradient | S11 | 154 | 0.403 | 0.936 | 0 | 0.936 / 0.921 |
| mid | S11 | 44 | 0.371 | 0.944 | 0 | 0.944 / 0.915 |
| held | S11v2 | 44 | 0.600 | 0.924 | 0 | 0.924 / 0.921 |
<!-- /table:fixed_hits -->

<!-- table:flips -->
| sample | run | boundary voxels with tries | voxels skipped | skipped voxels whose end the full pass had changed | kept voxels whose end differs from the full pass | fixed / broken vs the full pass | of those voxels, skipped ones | own hit ratio of those voxels |
|---|---|---|---|---|---|---|---|---|
| mid | S11n | 4068 | 13809 | 269 | 1674 | 4 / 7 | 0 | 0.885 to 0.915 |
| first gradient | S11n | 2968 | 20140 | 119 | 820 | 15 / 12 | 0 | 0.716 to 0.926 |
| big | S11 | 1828 | 9220 | 279 | 944 | 2 / 6 | 0 | 0.903 to 0.933 |
| held | S11v2 | 4497 | 16610 | 180 | 1264 | 4 / 10 | 0 | 0.770 to 0.915 |
<!-- /table:flips -->

<!-- table:live -->
| sample | run | pass evaluations, live | pass evaluations, replay | tries | voxels skipped | wrong, live (95% interval) | wrong, replay | wrong, full pass | fixed / broken vs the full pass | fixed / broken vs S8 | median error, live (deg) | run total, live | fewer than S8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid | S11n | 951,112 | 951,112 | 15592 (replay 15592) | 13809 (replay 13809) | 38 (0.11 to 0.21%) | 35 | 35 | 4 / 7 | 11 / 18 | 0.02913 | 11,997,038 | 54.6% |
| first gradient | S11n | 614,514 | 614,514 | 10074 (replay 10074) | 20140 (replay 20140) | 56 (0.18 to 0.30%) | 59 | 59 | 15 / 12 | 33 / 26 | 0.02123 | 6,469,435 | 64.1% |
| big | S11 | 273,951 | 273,951 | 4491 (replay 4491) | 9220 (replay 9220) | 22 (0.06 to 0.14%) | 18 | 18 | 2 / 6 | 9 / 7 | 0.03083 | 3,153,656 | 62.2% |
| held | S11v2 | 1,436,977 | 1,436,977 | 23557 (replay 23557) | 16610 (replay 16610) | 75 (0.24 to 0.38%) | 69 | 69 | 4 / 10 | 23 / 28 | 0.02616 | 33,753,877 | 43.7% |
<!-- /table:live -->

<!-- table:live_v2 -->
| sample | run | pass evaluations, live (equal to the replay) | kept tries identical to the stored full pass (voxel, candidate, end, cost, hit ratio, confidence, evaluations) | first difference | wrong, live (95% interval) | fixed / broken vs the full pass | fixed / broken vs S8 | median error (deg) | run total | fewer than S8 |
|---|---|---|---|---|---|---|---|---|---|---|
| mid | S11n | 951,112 | 15,592 of 15,592 | none | 35 (0.10 to 0.20%) | 0 / 0 | 13 / 17 | 0.02913 | 11,997,038 | 54.6% |
| first gradient | S11n | 614,514 | 10,074 of 10,074 | none | 59 (0.19 to 0.31%) | 0 / 0 | 33 / 29 | 0.02122 | 6,469,435 | 64.1% |
| big | S11 | 273,951 | 4,491 of 4,491 | none | 18 (0.05 to 0.12%) | 0 / 0 | 13 / 7 | 0.03084 | 3,153,656 | 62.2% |
| held (not clean: decided after held was read live under variant 1) | S11v2 | 1,436,977 | 23,557 of 23,557 | none | 69 (0.22 to 0.36%) | 0 / 0 | 24 / 23 | 0.02610 | 33,753,877 | 43.7% |
<!-- /table:live_v2 -->

<!-- table:noise -->
| sample | run | boundary voxels kept (with tries) | kept voxels whose end differs from the full pass | fixed / broken vs the full pass | wrong, live against full pass |
|---|---|---|---|---|---|
| mid | S11n | 4068 | 1674 | 4 / 7 | 38 against 35 |
| first gradient | S11n | 2968 | 820 | 15 / 12 | 56 against 59 |
| big | S11 | 1828 | 944 | 2 / 6 | 22 against 18 |
| held | S11v2 | 4497 | 1264 | 4 / 10 | 75 against 69 |
<!-- /table:noise -->

<!-- table:parity -->
| tries, live and stored | pass evaluations, live and stored | tries that differ: voxel / candidate / evaluations / cost / hit ratio / confidence / end | map equals the stored map | voxels changed, live and stored | same voxels changed |
|---|---|---|---|---|---|
| 50151 of 50151 | 3,059,211 of 3,059,211 | 0 / 0 / 0 / 0 / 0 / 0 / 0 | no: 1491 changed voxels differ by at most 5.7e-08 (matrix entries); 0 unchanged voxels differ | 1491 of 1491 | yes |
<!-- /table:parity -->

<!-- table:predictions -->
| prediction | sample | observed | verdict |
|---|---|---|---|
| 1 | held (replay) | cut 72.3%; fixed / broken 0 / 0 | met |
| 2 | held (live) | 10.0% fewer (33,753,877 against 37,507,939); the pass was 5,191,039 | met |
| 3 | mid S11n | live 951,112, replay 951,112; 15592 tries at 61.0 each | met |
| 3 | grad S11n | live 614,514, replay 614,514; 10074 tries at 61.0 each | met |
| 3 | big S11 | live 273,951, replay 273,951; 4491 tries at 61.0 each | met |
| 3 | held S11v2 | live 1,436,977, replay 1,436,977; 23557 tries at 61.0 each | met |
| 4 | mid S11n | live 38, replay 35 (difference +3) | not met |
| 4 | grad S11n | live 56, replay 59 (difference -3) | not met |
| 4 | big S11 | live 22, replay 18 (difference +4) | not met |
| 4 | held S11v2 | live 75, replay 69 (difference +6) | not met |
| 5 | mid S11n | +0.00000 degree against the full pass | met |
| 5 | grad S11n | +0.00005 degree against the full pass | met |
| 5 | big S11 | +0.00006 degree against the full pass | met |
| 5 | held S11v2 | +0.00011 degree against the full pass | met |
<!-- /table:predictions -->

<!-- table:predictions_v2 -->
| prediction | sample | observed | verdict |
|---|---|---|---|
| 1 | held (replay) | cut 72.3%; fixed / broken 0 / 0 | met |
| 2 | held (live) | 10.0% fewer (33,753,877 against 37,507,939); the pass was 5,191,039 | met (not clean: decided after held was read live under variant 1) |
| 3 | mid S11n | live 951,112, replay 951,112; 15592 tries at 61.0 each | met |
| 3 | grad S11n | live 614,514, replay 614,514; 10074 tries at 61.0 each | met |
| 3 | big S11 | live 273,951, replay 273,951; 4491 tries at 61.0 each | met |
| 3 | held S11v2 | live 1,436,977, replay 1,436,977; 23557 tries at 61.0 each | met (not clean: decided after held was read live under variant 1) |
| 4 | mid S11n | live 35, replay 35 (difference +0) | met |
| 4 | grad S11n | live 59, replay 59 (difference +0) | met |
| 4 | big S11 | live 18, replay 18 (difference +0) | met |
| 4 | held S11v2 | live 69, replay 69 (difference +0) | met (not clean: decided after held was read live under variant 1) |
| 5 | mid S11n | +0.00000 degree against the full pass | met |
| 5 | grad S11n | +0.00003 degree against the full pass | met |
| 5 | big S11 | +0.00007 degree against the full pass | met |
| 5 | held S11v2 | +0.00006 degree against the full pass | met (not clean: decided after held was read live under variant 1) |
<!-- /table:predictions_v2 -->

<!-- table:replay -->
| sample | run | disclosed | pass evaluations, full | pass evaluations, skip | cut | voxels skipped | wrong, full pass | wrong, skip (95% interval) | fixed / broken | median error, full (deg) | median error, skip (deg) | highest own hit ratio of a voxel wrong before the pass | rule parts of the pass met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid | S11n | yes | 3,059,211 | 951,112 | 68.9% | 13809 of 17877 | 35 | 35 (0.10 to 0.20%) | 0 / 0 | 0.02913 | 0.02913 | 0.944 | yes |
| first gradient | S11n | yes | 3,684,888 | 614,514 | 83.3% | 20140 of 23108 | 59 | 59 (0.19 to 0.31%) | 0 / 0 | 0.02118 | 0.02122 | 0.936 | yes |
| big | S11 | yes | 1,135,637 | 273,951 | 75.9% | 9220 of 11048 | 18 | 18 (0.05 to 0.12%) | 0 / 0 | 0.03077 | 0.03084 | 0.944 | yes |
| first gradient | S11 | yes | 3,685,803 | 609,634 | 83.5% | 20148 of 23101 | 62 | 62 (0.20 to 0.32%) | 0 / 0 | 0.02104 | 0.02105 | 0.936 | yes |
| mid | S11 | yes | 3,048,170 | 938,363 | 69.2% | 13821 of 17877 | 36 | 36 (0.11 to 0.20%) | 0 / 0 | 0.02913 | 0.02913 | 0.944 | yes |
| held | S11v2 | no (test) | 5,191,039 | 1,436,977 | 72.3% | 16610 of 21107 | 69 | 69 (0.22 to 0.36%) | 0 / 0 | 0.02604 | 0.02610 | 0.924 | yes |
<!-- /table:replay -->

<!-- table:rule_live -->
| sample | run | pass cut at least 50% | broken minus fixed vs the full pass within the allowance | median error at most 0.001 deg higher than the full pass | at least 20% fewer evaluations than S8 | broken minus fixed vs S8 at most 5 | median error at most 0.005 deg above S8 | the skip passes |
|---|---|---|---|---|---|---|---|---|
| mid | S11n | met (68.9%) | not met (3 <= 2) | met (+0.00000) | met | not met (7 <= 5) | met (-0.00060) | no |
| first gradient | S11n | met (83.3%) | met (-3 <= 8) | met (+0.00005) | met | met (-7 <= 5) | met (-0.00141) | yes |
| big | S11 | met (75.9%) | not met (4 <= 3) | met (+0.00006) | met | met (-2 <= 5) | met (-0.00078) | no |
| held | S11v2 | met (72.3%) | not met (6 <= 2) | met (+0.00011) | met | met (5 <= 5) | met (-0.00243) | no |
<!-- /table:rule_live -->

<!-- table:rule_live_v2 -->
| sample | run | pass cut at least 50% | broken minus fixed vs the full pass within the allowance | median error at most 0.001 deg higher than the full pass | at least 20% fewer evaluations than S8 | broken minus fixed vs S8 at most 5 | median error at most 0.005 deg above S8 | the skip passes |
|---|---|---|---|---|---|---|---|---|
| mid | S11n | met (68.9%) | met (0 <= 2) | met (+0.00000) | met | met (4 <= 5) | met (-0.00060) | yes |
| first gradient | S11n | met (83.3%) | met (0 <= 8) | met (+0.00003) | met | met (-4 <= 5) | met (-0.00142) | yes |
| big | S11 | met (75.9%) | met (0 <= 3) | met (+0.00007) | met | met (-6 <= 5) | met (-0.00077) | yes |
| held | S11v2 | met (72.3%) | met (0 <= 2) | met (+0.00006) | met | met (-1 <= 5) | met (-0.00248) | yes (not clean: decided after held was read live under variant 1) |
<!-- /table:rule_live_v2 -->

<!-- table:rule_replay -->
| sample | run | pass cut at least 50% | broken minus fixed within the allowance | median within 0.001 deg | house rules vs S8 | the skip passes (replay) |
|---|---|---|---|---|---|---|
| mid | S11n | met | met | met | met | yes |
| first gradient | S11n | met | met | met | met | yes |
| big | S11 | met | met | met | met | yes |
| first gradient | S11 | met | met | met | met | yes |
| mid | S11 | met | met | met | met | yes |
| held | S11v2 | met | met | met | met | yes |
<!-- /table:rule_replay -->

<!-- table:variants -->
| sample | run | variant 1 (one draw stream): passes the rule | variant 1: fixed / broken vs the full pass | variant 2 (aligned draws): passes the rule | variant 2: fixed / broken vs the full pass |
|---|---|---|---|---|---|
| mid | S11n | no | 4 / 7 | yes | 0 / 0 |
| first gradient | S11n | yes | 15 / 12 | yes | 0 / 0 |
| big | S11 | no | 2 / 6 | yes | 0 / 0 |
| held | S11v2 | no | 4 / 10 | yes (not clean: decided after held was read live under variant 1) | 0 / 0 |
<!-- /table:variants -->
