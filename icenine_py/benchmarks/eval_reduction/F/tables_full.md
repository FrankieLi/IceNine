<!-- table:evalF_full_acc -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F4 Gauss-Newton on E4 rays | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 81 accepted, 0 of them wrong; 0 of 39 rejected were right |
| F5 expanding shells | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right |
| F5 then F4 | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right |
<!-- /table:evalF_full_acc -->

<!-- table:evalF_full_acc_margin -->
| method | runs | accepted but wrong | rejected but right | largest accepted strict cost | smallest rejected strict cost | largest error accepted (deg) | smallest error rejected (deg) |
|---|---|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | all offsets | 0 of 321 (0.0-1.2%) | 0 of 39 (0.0-9.0%) | 0.437 | 0.953 | 0.295 | 3.37 |
| F4 Gauss-Newton on E4 rays | 5 deg | 0 of 81 (0.0-4.5%) | 0 of 39 (0.0-9.0%) | 0.159 | 0.953 | 0.086 | 3.37 |
| F4 with a 5-frame matching window | all offsets | 0 of 329 (0.0-1.2%) | 0 of 31 (0.0-11.0%) | 0.437 | 0.953 | 0.295 | 3.42 |
| F4 with a 5-frame matching window | 5 deg | 0 of 89 (0.0-4.1%) | 0 of 31 (0.0-11.0%) | 0.264 | 0.953 | 0.244 | 3.42 |
<!-- /table:evalF_full_acc_margin -->

<!-- table:evalF_full_cap005 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 81% (73-87) | 63% (54-71) | 22% (15-30) |
| F2 CMA-ES, sigma0 0.5 deg | 89% (82-94) | 82% (74-88) | 8% (5-15) |
| F2 CMA-ES, sigma0 1 deg | 92% (85-95) | 89% (82-94) | 41% (32-50) |
| F2 CMA-ES, sigma0 2 deg | 91% (84-95) | 92% (86-96) | 66% (57-74) |
| F4 Gauss-Newton on E4 rays | 33% (26-42) | 33% (26-42) | 25% (18-33) |
| F4 with a 5-frame matching window | 50% (41-59) | 50% (41-59) | 37% (29-46) |
| F5 expanding shells | 0% (0-3) | 0% (0-3) | 0% (0-3) |
| F5 then F4 | 33% (26-42) | 33% (26-42) | 33% (26-42) |
<!-- /table:evalF_full_cap005 -->

<!-- table:evalF_full_cap01 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 100% (97-100) | 87% (79-92) | 32% (25-41) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (97-100) | 92% (86-96) | 10% (6-17) |
| F2 CMA-ES, sigma0 1 deg | 100% (97-100) | 100% (97-100) | 42% (33-51) |
| F2 CMA-ES, sigma0 2 deg | 97% (92-99) | 97% (92-99) | 73% (65-80) |
| F4 Gauss-Newton on E4 rays | 98% (94-100) | 99% (95-100) | 68% (59-75) |
| F4 with a 5-frame matching window | 98% (94-100) | 100% (97-100) | 72% (63-79) |
| F5 expanding shells | 0% (0-3) | 0% (0-3) | 0% (0-3) |
| F5 then F4 | 98% (94-100) | 98% (94-100) | 99% (95-100) |
<!-- /table:evalF_full_cap01 -->

<!-- table:evalF_full_cap1 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 100% (97-100) | 96% (91-98) | 56% (47-64) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (97-100) | 92% (86-96) | 11% (6-18) |
| F2 CMA-ES, sigma0 1 deg | 100% (97-100) | 100% (97-100) | 42% (33-51) |
| F2 CMA-ES, sigma0 2 deg | 97% (92-99) | 98% (93-99) | 76% (67-83) |
| F4 Gauss-Newton on E4 rays | 100% (97-100) | 100% (97-100) | 68% (59-75) |
| F4 with a 5-frame matching window | 100% (97-100) | 100% (97-100) | 74% (66-81) |
| F5 expanding shells | 100% (97-100) | 100% (97-100) | 100% (97-100) |
| F5 then F4 | 100% (97-100) | 100% (97-100) | 100% (97-100) |
<!-- /table:evalF_full_cap1 -->

<!-- table:evalF_full_choice -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 2.99 | 1.82 | 3 | 2.4 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.2 | 0.95 | 4 | 2.08 | yes |
| F4 Gauss-Newton on E4 rays | 2.69 | 1.52 | 6 | 2.1 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.86 | 0.276 | 10 | 1.57 |  |
| F5 expanding shells | 2.29 | 0.494 | 10 | 1.39 |  |
| F1 FindOptimal, box = offset | 2.05 | 0.642 | 10 | 1.35 |  |
| F5 then F4 | 1.41 | 0.436 | 13 | 0.923 |  |
<!-- /table:evalF_full_choice -->

<!-- table:evalF_full_choice_v1 -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 2.99 | 1.82 | 3 | 2.4 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.2 | 0.95 | 4 | 2.08 | yes |
| F4 Gauss-Newton on E4 rays | 2.69 | 1.52 | 6 | 2.1 |  |
| F1 FindOptimal, box = offset | 2.05 | 0.642 | 9 | 1.35 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.86 | 0.276 | 10 | 1.57 |  |
| F5 expanding shells | 1.86 | 0.396 | 11 | 1.13 |  |
| F5 then F4 | 1.24 | 0.359 | 13 | 0.797 |  |
<!-- /table:evalF_full_choice_v1 -->

<!-- table:evalF_full_choice_wall -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 3.04 | 1.83 | 3 | 2.44 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.23 | 0.966 | 4 | 2.1 | yes |
| F4 Gauss-Newton on E4 rays | 2.86 | 1.6 | 6 | 2.23 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.87 | 0.28 | 10 | 1.58 |  |
| F1 FindOptimal, box = offset | 2.14 | 0.669 | 10 | 1.4 |  |
| F5 expanding shells | 2.21 | 0.484 | 10 | 1.35 |  |
| F5 then F4 | 1.4 | 0.431 | 13 | 0.917 |  |
<!-- /table:evalF_full_choice_wall -->

<!-- table:evalF_full_choice_with_f4w -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 2.99 | 1.82 | 3 | 2.4 |  |
| F2 CMA-ES, sigma0 1 deg | 3.2 | 0.95 | 5 | 2.08 |  |
| F4 with a 5-frame matching window | 2.72 | 1.61 | 6 | 2.17 |  |
| F4 Gauss-Newton on E4 rays | 2.69 | 1.52 | 8 | 2.1 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.86 | 0.276 | 11 | 1.57 |  |
| F5 expanding shells | 2.29 | 0.494 | 12 | 1.39 |  |
| F1 FindOptimal, box = offset | 2.05 | 0.642 | 12 | 1.35 |  |
| F5 then F4 | 1.41 | 0.436 | 15 | 0.923 |  |
<!-- /table:evalF_full_choice_with_f4w -->

<!-- table:evalF_full_err -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 0.0297 | 0.0371 | 0.4314 |
| F2 CMA-ES, sigma0 0.5 deg | 0.0274 | 0.0275 | 4.8527 |
| F2 CMA-ES, sigma0 1 deg | 0.0242 | 0.0276 | 3.8824 |
| F2 CMA-ES, sigma0 2 deg | 0.0255 | 0.0278 | 0.0353 |
| F4 Gauss-Newton on E4 rays | 0.0520 | 0.0520 | 0.0538 |
| F4 with a 5-frame matching window | 0.0486 | 0.0486 | 0.0519 |
| F5 expanding shells | 0.4391 | 0.5281 | 0.5176 |
| F5 then F4 | 0.0520 | 0.0514 | 0.0519 |
<!-- /table:evalF_full_err -->

<!-- table:evalF_full_pervoxel_cap1_5 -->
| method | voxel 5242 | voxel 15901 | voxel 17034 | voxel 18558 | voxel 19889 | voxel 22867 |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 11 of 20 | 12 of 20 | 10 of 20 | 15 of 20 | 7 of 20 | 12 of 20 |
| F2 CMA-ES, sigma0 1 deg | 9 of 20 | 9 of 20 | 6 of 20 | 7 of 20 | 7 of 20 | 12 of 20 |
| F2 CMA-ES, sigma0 2 deg | 15 of 20 | 14 of 20 | 18 of 20 | 15 of 20 | 15 of 20 | 14 of 20 |
| F4 Gauss-Newton on E4 rays | 14 of 20 | 18 of 20 | 14 of 20 | 12 of 20 | 13 of 20 | 10 of 20 |
| F4 with a 5-frame matching window | 14 of 20 | 19 of 20 | 17 of 20 | 14 of 20 | 14 of 20 | 11 of 20 |
| F5 expanding shells | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 |
| F5 then F4 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 |
<!-- /table:evalF_full_pervoxel_cap1_5 -->

<!-- table:evalF_full_policy -->
| offset | runs | F4 rejected, F5+F4 run | within 1 deg | within 0.1 deg | within 0.05 deg | median evals (strict-equivalent) | 90th percentile |
|---|---|---|---|---|---|---|---|
| 1 deg | 120 | 0 | 100% (97-100) | 98% (94-100) | 33% (26-42) | 317 | 431 |
| 2 deg | 120 | 0 | 100% (97-100) | 99% (95-100) | 33% (26-42) | 372 | 499 |
| 5 deg | 120 | 39 | 100% (97-100) | 99% (95-100) | 33% (26-42) | 547 | 2806 |
<!-- /table:evalF_full_policy -->

<!-- table:evalF_full_policy_f4w -->
| offset | runs | F4 w5 rejected, F5+F4 run | within 1 deg | within 0.1 deg | within 0.05 deg | median evals (strict-equivalent) | 90th percentile |
|---|---|---|---|---|---|---|---|
| 1 deg | 120 | 0 | 100% (97-100) | 98% (94-100) | 50% (41-59) | 319 | 403 |
| 2 deg | 120 | 0 | 100% (97-100) | 100% (97-100) | 50% (41-59) | 367 | 465 |
| 5 deg | 120 | 31 | 100% (97-100) | 98% (93-99) | 42% (34-51) | 502 | 2816 |
<!-- /table:evalF_full_policy_f4w -->

<!-- table:evalF_full_se_vs_wall -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 1.05 -> 1.05 | 1.04 -> 1.04 | 1.04 -> 1.04 |
| F2 CMA-ES, sigma0 0.5 deg | 1.02 -> 1.02 | 1.00 -> 1.00 | 1.01 -> 1.01 |
| F2 CMA-ES, sigma0 1 deg | 1.01 -> 1.01 | 1.01 -> 1.01 | 1.02 -> 1.02 |
| F2 CMA-ES, sigma0 2 deg | 1.02 -> 1.02 | 1.02 -> 1.02 | 1.01 -> 1.01 |
| F4 Gauss-Newton on E4 rays | 1.07 -> 1.07 | 1.06 -> 1.06 | 1.05 -> 1.05 |
| F4 with a 5-frame matching window | 1.06 -> 1.06 | 1.06 -> 1.06 | 1.05 -> 1.05 |
| F5 expanding shells | 1.18 -> 0.96 | 1.19 -> 0.97 | 1.22 -> 0.98 |
| F5 then F4 | 1.10 -> 1.02 | 1.15 -> 0.99 | 1.20 -> 0.99 |
<!-- /table:evalF_full_se_vs_wall -->

<!-- table:evalF_full_secs -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 0.223 | 0.219 | 0.390 |
| F2 CMA-ES, sigma0 0.5 deg | 0.134 | 0.153 | 0.187 |
| F2 CMA-ES, sigma0 1 deg | 0.132 | 0.146 | 0.202 |
| F2 CMA-ES, sigma0 2 deg | 0.157 | 0.154 | 0.201 |
| F4 Gauss-Newton on E4 rays | 0.142 | 0.160 | 0.200 |
| F4 with a 5-frame matching window | 0.147 | 0.160 | 0.200 |
| F5 expanding shells | 0.065 | 0.206 | 0.949 |
| F5 then F4 | 0.181 | 0.325 | 1.072 |
<!-- /table:evalF_full_secs -->

<!-- table:evalF_full_secs_single -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 0.202 | 0.215 | 0.367 |
| F2 CMA-ES, sigma0 0.5 deg | 0.138 | 0.145 | 0.146 |
| F2 CMA-ES, sigma0 1 deg | 0.132 | 0.152 | 0.207 |
| F2 CMA-ES, sigma0 2 deg | 0.155 | 0.154 | 0.221 |
| F4 Gauss-Newton on E4 rays | 0.143 | 0.150 | 0.200 |
| F4 with a 5-frame matching window | 0.132 | 0.144 | 0.185 |
| F5 expanding shells | 0.062 | 0.199 | 0.909 |
| F5 then F4 | 0.160 | 0.293 | 1.013 |
<!-- /table:evalF_full_secs_single -->

<!-- table:evalF_full_stop_raw -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 496 (540) | 467 (553) | 870 (902) |
| F2 CMA-ES, sigma0 0.5 deg | 281 (422) | 323 (527) | 393 (569) |
| F2 CMA-ES, sigma0 1 deg | 288 (415) | 312 (478) | 438 (625) |
| F2 CMA-ES, sigma0 2 deg | 334 (499) | 326 (464) | 418 (682) |
| F4 Gauss-Newton on E4 rays | 29 (39.1) | 31.5 (43) | 38 (57) |
| F4 with a 5-frame matching window | 29 (36) | 31 (39.1) | 38 (53.1) |
| F5 expanding shells | 125 (125) | 954 (955) | 18495 (18496) |
| F5 then F4 | 145 (156) | 979 (1001) | 18520 (18539) |
<!-- /table:evalF_full_stop_raw -->

<!-- table:evalF_full_stop_se -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 496 (540) | 467 (553) | 870 (902) |
| F2 CMA-ES, sigma0 0.5 deg | 281 (422) | 323 (527) | 393 (569) |
| F2 CMA-ES, sigma0 1 deg | 288 (415) | 312 (478) | 438 (625) |
| F2 CMA-ES, sigma0 2 deg | 334 (499) | 326 (464) | 418 (682) |
| F4 Gauss-Newton on E4 rays | 317 (431) | 372 (499) | 444 (658) |
| F4 with a 5-frame matching window | 319 (403) | 367 (465) | 461 (657) |
| F5 expanding shells | 137 (137) | 437 (438) | 2026 (2027) |
| F5 then F4 | 392 (474) | 709 (813) | 2292 (2436) |
<!-- /table:evalF_full_stop_se -->

<!-- table:evalF_full_timing -->
| kind of evaluation | used for | milliseconds | times a strict evaluation |
|---|---|---|---|
| strict cost | all methods | 0.4398 | 1.000 |
| screen (pixel radius 3), per orientation, batches of 64 | v1 weights only | 0.0224 | 0.051 |
| E4 value, one rotation per call | F4 | 1.4978 | 3.406 |
| E4 value and Jacobian, one rotation per call | F4 | 12.0732 | 27.452 |
| screen, per orientation, a shell at a time (chunks of up to 256) | F5 | 0.0175 | 0.040 |
| E4 value, per orientation, 32 rotations per call | F5 | 1.2144 | 2.761 |
<!-- /table:evalF_full_timing -->

<!-- table:evalF_full_timing_voxels -->
| voxel | weights from it | F4 E4 value (ms) | F4 E4 value and Jacobian (ms) | F5 screen (ms per orientation) | F5 E4 (ms per orientation) | strict inside F5 (ms) |
|---|---|---|---|---|---|---|
| 15901 | yes | 1.456 | 11.63 | 0.0175 | 1.214 | 0.506 |
<!-- /table:evalF_full_timing_voxels -->

<!-- table:evalF_full_tocap005 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 75 (102) n=107 | 112 (184) n=80 | 168 (251) n=28 |
| F2 CMA-ES, sigma0 0.5 deg | 105 (151) n=114 | 146 (194) n=103 | 265 (352) n=11 |
| F2 CMA-ES, sigma0 1 deg | 118 (156) n=114 | 147 (190) n=113 | 239 (446) n=49 |
| F2 CMA-ES, sigma0 2 deg | 153 (208) n=111 | 158 (220) n=111 | 220 (296) n=82 |
| F4 Gauss-Newton on E4 rays | 15 (21) n=73 | 19 (25) n=77 | 25 (36.4) n=54 |
| F4 with a 5-frame matching window | 15 (25.4) n=84 | 19 (28) n=90 | 27 (37) n=66 |
| F5 expanding shells | n/a | n/a | n/a |
| F5 then F4 | 136 (141) n=70 | 967 (1518) n=70 | 18508 (18516) n=75 |
<!-- /table:evalF_full_tocap005 -->

<!-- table:evalF_full_tocap1 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 5 (15.1) n=120 | 23 (166) n=115 | 107 (244) n=67 |
| F2 CMA-ES, sigma0 0.5 deg | 4 (11) n=120 | 26 (73) n=111 | 111 (206) n=13 |
| F2 CMA-ES, sigma0 1 deg | 10 (29) n=120 | 18 (54) n=120 | 95 (256) n=50 |
| F2 CMA-ES, sigma0 2 deg | 26 (64.5) n=116 | 27 (65.8) n=117 | 81 (176) n=91 |
| F4 Gauss-Newton on E4 rays | 5 (5) n=120 | 7 (9) n=120 | 13 (21) n=81 |
| F4 with a 5-frame matching window | 5 (5) n=120 | 7 (7) n=120 | 13 (21) n=89 |
| F5 expanding shells | 33 (33) n=120 | 558 (558) n=120 | 15197 (15197) n=120 |
| F5 then F4 | 33 (33) n=120 | 558 (558) n=120 | 15197 (15197) n=120 |
<!-- /table:evalF_full_tocap1 -->
