<!-- table:ar1_broken -->
| voxel | taken | boundary | cost_taken | cost_truth | fin_taken | fin_truth | start_deg | against_best1 |
|---|---|---|---|---|---|---|---|---|
| 7141 | grain 23 | yes (shares an edge) | 0.183 | 0.201 | 0.177 | 0.080 | 49.5 | yes |
| 8848 | grain 395 | yes (shares an edge) | 0.187 | 0.190 | 0.139 | 0.028 | 43.8 | yes |
| 23942 | grain 286 | yes (shares an edge) | 0.271 | 0.335 | 0.209 | 0.111 | 54.4 | yes |
<!-- /table:ar1_broken -->

<!-- table:ar1_dist -->
| set | best_med | best_p90 | best_max | top1_med | n_cand_med |
|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 0.038 | 0.076 | 0.125 | 0.038 | 20 |
| Pilot full searches (69 voxels) | 0.054 | 0.082 | 0.137 | 0.048 | 20 |
<!-- /table:ar1_dist -->

<!-- table:ar1_e2e -->
| set | arm | wrong | fixed | broken | both | mean_ev | med_ev | med_err |
|---|---|---|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | stored full search | 41/200 (20.5%, 15.5-26.6) | 0 | 0 | 41 | 48,700 | 48,724 | 0.021 |
| E0 (200 voxels, easy: own spots only) | strict cost on first 1, no finisher | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 1 | 1 | 0.038 |
| E0 (200 voxels, easy: own spots only) | strict cost on first 5, no finisher | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 5 | 5 | 0.038 |
| E0 (200 voxels, easy: own spots only) | strict cost on first 20, no finisher | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 20 | 20 | 0.038 |
| E0 (200 voxels, easy: own spots only) | AR1: strict on 5, finisher on best 1 (set in advance) | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 519 | 530 | 0.017 |
| E0 (200 voxels, easy: own spots only) | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 950 | 951 | 0.017 |
| E0 (200 voxels, easy: own spots only) | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 1,393 | 1,392 | 0.017 |
| Pilot full searches (69 voxels) | stored full search | 4/69 (5.8%, 2.3-14.0) | 0 | 0 | 4 | 332,380 | 332,154 | 0.015 |
| Pilot full searches (69 voxels) | strict cost on first 1, no finisher | 28/69 (40.6%, 29.8-52.4) | 2 | 26 | 2 | 1 | 1 | 0.073 |
| Pilot full searches (69 voxels) | strict cost on first 5, no finisher | 3/69 (4.3%, 1.5-12.0) | 4 | 3 | 0 | 5 | 5 | 0.054 |
| Pilot full searches (69 voxels) | strict cost on first 20, no finisher | 3/69 (4.3%, 1.5-12.0) | 4 | 3 | 0 | 20 | 20 | 0.054 |
| Pilot full searches (69 voxels) | AR1: strict on 5, finisher on best 1 (set in advance) | 3/69 (4.3%, 1.5-12.0) | 4 | 3 | 0 | 1,125 | 1,138 | 0.020 |
| Pilot full searches (69 voxels) | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/69 (0.0%, 0.0-5.3) | 4 | 0 | 0 | 1,926 | 1,907 | 0.019 |
| Pilot full searches (69 voxels) | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/69 (0.0%, 0.0-5.3) | 4 | 0 | 0 | 2,546 | 2,512 | 0.019 |
<!-- /table:ar1_e2e -->

<!-- table:ar1_e2e_split -->
| subset | arm | wrong | fixed | broken | both |
|---|---|---|---|---|---|
| 62 seeds | stored full search | 2/62 (3.2%, 0.9-11.0) | 0 | 0 | 2 |
| 62 seeds | AR1: strict on 5, finisher on best 1 (set in advance) | 3/62 (4.8%, 1.7-13.3) | 2 | 3 | 0 |
| 62 seeds | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/62 (0.0%, 0.0-5.8) | 2 | 0 | 0 |
| 62 seeds | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/62 (0.0%, 0.0-5.8) | 2 | 0 | 0 |
| 7 rejected seeds | stored full search | 2/7 (28.6%, 8.2-64.1) | 0 | 0 | 2 |
| 7 rejected seeds | AR1: strict on 5, finisher on best 1 (set in advance) | 0/7 (0.0%, 0.0-35.4) | 2 | 0 | 0 |
| 7 rejected seeds | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/7 (0.0%, 0.0-35.4) | 2 | 0 | 0 |
| 7 rejected seeds | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/7 (0.0%, 0.0-35.4) | 2 | 0 | 0 |
<!-- /table:ar1_e2e_split -->

<!-- table:ar1_index_time -->
| set | med_ms | p90_ms | med_votes | med_att | med_true |
|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 91 | 101 | 5022 | 57 | 57 |
| Pilot full searches (69 voxels) | 179 | 234 | 28248 | 413 | 44 |
<!-- /table:ar1_index_time -->

<!-- table:ar1_misses -->
| what | n |
|---|---|
| rank-1 candidate is not the truth | 28/69 |
|   reason: too few rays (truth explains < 8) | 0 |
|   reason: another real grain explains more rays | 28 |
|   reason: matches no grain (coincidence of rays) | 0 |
| rank-1 grain touches the voxel (shares a vertex / an edge) | 27 / 16 |
| truth is candidate 1 / 2 / 3 (of 20) | 41 / 22 / 6 |
<!-- /table:ar1_misses -->

<!-- table:ar1_ray_kinds -->
| quantity | median |
|---|---|
| attached rays | 413 |
| drawn only by the voxel's own grain | 45 |
| drawn by the own grain and another grain | 1 |
| drawn by other grains only | 366 |
| explained by the true orientation | 44 |
<!-- /table:ar1_ray_kinds -->

<!-- table:ar1_recall -->
| set | K1 | K5 | K20 |
|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 200/200 (100.0%, 98.1-100.0) | 200/200 (100.0%, 98.1-100.0) | 200/200 (100.0%, 98.1-100.0) |
| Pilot full searches (69 voxels) | 41/69 (59.4%, 47.6-70.2) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
<!-- /table:ar1_recall -->

<!-- table:ar1_recall_thr -->
| set | within | K1 | K5 | K20 |
|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 0.5 deg | 200/200 | 200/200 | 200/200 |
| E0 (200 voxels, easy: own spots only) | 2 deg | 200/200 | 200/200 | 200/200 |
| E0 (200 voxels, easy: own spots only) | 5 deg | 200/200 | 200/200 | 200/200 |
| Pilot full searches (69 voxels) | 0.5 deg | 41/69 | 69/69 | 69/69 |
| Pilot full searches (69 voxels) | 2 deg | 41/69 | 69/69 | 69/69 |
| Pilot full searches (69 voxels) | 5 deg | 41/69 | 69/69 | 69/69 |
<!-- /table:ar1_recall_thr -->

<!-- table:ar1_sanity -->
| set | voxels | attached | t05 | c05 | t1 | c1 | ang50 | ang95 | angmax | len95 |
|---|---|---|---|---|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 200 | 11312 | 11270 (99.6%) | 3 (0.03%) | 11312 | 14 | 0.203 | 0.439 | 0.559 | 0.0026 |
| Pilot full searches (69 voxels) | 69 | 27164 | 3065 (11.3%) | 4 (0.01%) | 3115 | 28 | 0.206 | 0.451 | 0.985 | 0.0056 |
<!-- /table:ar1_sanity -->

<!-- table:ar1_sparse -->
| case | rays | K1 | K5 | K20 |
|---|---|---|---|---|
| stride 1 | 44 | 41/69 (59.4%, 47.6-70.2) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| stride 2 | 23 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| stride 3 | 15 | 37/69 (53.6%, 42.0-64.9) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| stride 4 | 11 | 33/69 (47.8%, 36.5-59.4) | 67/69 (97.1%, 90.0-99.2) | 68/69 (98.6%, 92.2-99.7) |
| stride 6 | 7 | 32/69 (46.4%, 35.1-58.0) | 62/69 (89.9%, 80.5-95.0) | 67/69 (97.1%, 90.0-99.2) |
| stride 8 | 5 | 30/69 (43.5%, 32.4-55.2) | 54/69 (78.3%, 67.2-86.4) | 66/69 (95.7%, 88.0-98.5) |
| ang 0.35 | 36 | 44/69 (63.8%, 52.0-74.1) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| ang 1.4 | 46 | 41/69 (59.4%, 47.6-70.2) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| pair 0.5 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| pair 2 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| cluster 0.5 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| cluster 2 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
<!-- /table:ar1_sparse -->

<!-- table:ar1_timing -->
| quantity | seconds |
|---|---|
| one-time pairing of all spots (26,019 rays) | 17.069 |
| attach rays to a voxel (median) | 0.002 |
| index a voxel (median) | 0.177 |
| index a voxel (90th percentile) | 0.230 |
|   of which voting (median) | 0.021 |
|   of which clustering (median) | 0.085 |
|   of which refining (median) | 0.069 |
<!-- /table:ar1_timing -->
