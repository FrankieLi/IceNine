<!-- table:gap -->
| aspect | measured | deployed | gap |
|---|---|---|---|
| success criterion | final error < 1 deg from a 1/2/5 deg start; recorded 96/41/0% (Adam), 92/40/6% (MC). At a 1 deg start this means 'ended below the start' | a finisher meant to end within about 0.01-0.03 deg of the truth | the criterion does not measure the deployed target |
| precision of the 'successes' | median final error 0.47 deg (Adam) and 0.53 deg (MC); under 0.1 deg: 3/100 and 5/100; under 0.02 deg: 0/100 and 0/100 | finisher result median 0.0229 deg from the truth (T5, H3 realistic) | the sweep's successes are about 23x (MC) coarser than the finisher's result; the 0.01-0.03 deg scale was not measured |
| MC settings | 3500 steps, 2 restarts, step frac 0.5; box = 1.5 x the start offset, step = 0.5 x box (0.75 deg at a 1 deg start) | 200 steps, 2 restarts, box 0.329 deg, step 0.1317 deg | different step scale, box and budget |
| starting point | the truth rotated by exactly 1, 2 or 5 deg about a random axis | a coarse-search or network start; T5 H3 start error 0.088 deg (q25-q75 0.030-0.152) | the sweep's starts are 11x (1 deg) to 57x (5 deg) the T5 H3 median start error |
| budget | Adam/SGD 101-501 evaluations (gradient evaluations); MC 101-3501; hybrid 4 hard + 303 differentiable | MC 200 steps, then VarianceMinimizing; whole finisher median 2629 evaluations (T5 H3) | MC at 101 and at 3501 evaluations give the same result (next row) |
| matched budget (MC) | 100 steps (best config): 49/100 (49%, 39-59) under 0.5 deg at a 1 deg start, median error of the < 1 deg runs 0.465; 3500 steps (best config, 5 restarts): 50/100 (50%, 40-60), 0.461; pooled over each class 242/900 (27%, 24-30) vs 245/900 (27%, 24-30); the recorded 3500-step 2-restart config: 39/100 (39%, 30-49), 0.531 | 200 steps (not run in the sweep) | more MC steps bought no accuracy in the sweep between 100 and 3500 steps; 200 steps were not run |
| data and voxels | ManyGrains: 100 voxels, r_perp 75-563 um (median 370); ThreeVoxels: 3 voxels at 12 um | per-voxel finisher on the reconstruction's voxels | no detectable r_perp dependence inside the sweep's range (Spearman of error vs r_perp, recorded Adam config: rho 0.02, p 0.85, n 100 voxels) |
| MC stopping in the sweep | recorded MC config at a 1 deg start: last accepted move before 10% of the budget in 92/99 (93%, 86-97) | same rule, 200 steps | consistent with the step collapse measured in B1; never examined in the sweep |
<!-- /table:gap -->

<!-- table:hybrid -->
| example | start offset (deg) | n | hybrid < 0.5 deg | MC < 0.5 deg | hybrid < 1 deg | MC < 1 deg | hybrid median error of < 0.5 deg runs | MC median error of < 0.5 deg runs | hybrid median error, all | MC median error, all | hybrid evaluations | MC evaluations | hybrid time (s) | MC time (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| manygrains | 0.5 | 20 | 20/20 (100%, 84-100) | 18/20 (90%, 70-97) | 20/20 (100%, 84-100) | 20/20 (100%, 84-100) | 0.224 | 0.249 | 0.224 | 0.259 | 4+303d | 3501 | 0.81 | 1.88 |
| manygrains | 1 | 20 | 12/20 (60%, 39-78) | 9/20 (45%, 26-66) | 20/20 (100%, 84-100) | 20/20 (100%, 84-100) | 0.346 | 0.295 | 0.408 | 0.592 | 4+303d | 3501 | 0.81 | 1.88 |
| manygrains | 2 | 20 | 1/20 (5%, 1-24) | 3/20 (15%, 5-36) | 5/20 (25%, 11-47) | 7/20 (35%, 18-57) | 0.355 | 0.078 | 1.101 | 1.429 | 4+303d | 3501 | 0.81 | 1.88 |
| manygrains | 3 | 20 | 0/20 (0%, 0-16) | 2/20 (10%, 3-30) | 0/20 (0%, 0-16) | 4/20 (20%, 8-42) | n/a | 0.354 | 3.000 | 3.093 | 4+303d | 3501 | 0.81 | 1.88 |
| manygrains | 5 | 20 | 0/20 (0%, 0-16) | 0/20 (0%, 0-16) | 0/20 (0%, 0-16) | 0/20 (0%, 0-16) | n/a | n/a | 5.431 | 7.178 | 4+303d | 3501 | 0.81 | 1.86 |
| threevoxels | 0.5 | 3 | 3/3 (100%, 44-100) | 2/3 (67%, 21-94) | 3/3 (100%, 44-100) | 3/3 (100%, 44-100) | 0.274 | 0.229 | 0.274 | 0.259 | 4+303d | 3501 | 3.15 | 4.77 |
| threevoxels | 1 | 3 | 2/3 (67%, 21-94) | 2/3 (67%, 21-94) | 3/3 (100%, 44-100) | 3/3 (100%, 44-100) | 0.306 | 0.376 | 0.463 | 0.459 | 4+303d | 3501 | 2.97 | 4.93 |
| threevoxels | 2 | 3 | 0/3 (0%, 0-56) | 0/3 (0%, 0-56) | 0/3 (0%, 0-56) | 1/3 (33%, 6-79) | n/a | n/a | 2.000 | 2.000 | 4+303d | 3501 | 2.93 | 4.79 |
| threevoxels | 3 | 3 | 0/3 (0%, 0-56) | 0/3 (0%, 0-56) | 1/3 (33%, 6-79) | 0/3 (0%, 0-56) | n/a | n/a | 2.638 | 3.170 | 4+303d | 3501 | 2.88 | 4.87 |
| threevoxels | 5 | 3 | 0/3 (0%, 0-56) | 0/3 (0%, 0-56) | 0/3 (0%, 0-56) | 0/3 (0%, 0-56) | n/a | n/a | 5.000 | 4.635 | 4+303d | 3501 | 2.93 | 4.78 |
<!-- /table:hybrid -->

<!-- table:hybrid_precision -->
| example | optimizer | runs (offsets <= 1 deg) | final error q25/50/75 | minimum | < 0.5 deg | < 0.1 deg | < 0.02 deg |
|---|---|---|---|---|---|---|---|
| manygrains | hybrid | 40 | 0.168/0.326/0.436 | 0.079 | 32/40 (80%, 65-90) | 2/40 (5%, 1-17) | 0/40 (0%, 0-9) |
| manygrains | mc | 40 | 0.237/0.351/0.589 | 0.057 | 27/40 (68%, 52-80) | 2/40 (5%, 1-17) | 0/40 (0%, 0-9) |
| threevoxels | hybrid | 6 | 0.180/0.350/0.454 | 0.028 | 5/6 (83%, 44-97) | 1/6 (17%, 3-56) | 0/6 (0%, 0-39) |
| threevoxels | mc | 6 | 0.267/0.376/0.586 | 0.199 | 4/6 (67%, 30-90) | 0/6 (0%, 0-39) | 0/6 (0%, 0-39) |
<!-- /table:hybrid_precision -->

<!-- table:mc_headline_timing -->
| start offset (deg) | n | ended before the budget | median accepts | step of last accept q25/50/75 | last accept before 10% of budget | step at last accept (deg) q25/50/75 | median final error (deg) |
|---|---|---|---|---|---|---|---|
| 1 | 100 | 2/100 (2%, 1-7) | 10 | 27/40/70 | 92/99 (93%, 86-97) | 0.0004/0.0007/0.0015 | 0.579 |
| 2 | 100 | 2/100 (2%, 1-7) | 9 | 30/42/72 | 95/98 (97%, 91-99) | 0.0007/0.0029/0.0059 | 1.339 |
| 5 | 100 | 5/100 (5%, 2-11) | 9 | 62/115/262 | 81/97 (84%, 75-90) | 0.0018/0.0073/0.0146 | 5.368 |
<!-- /table:mc_headline_timing -->

<!-- table:mc_stopping -->
| max steps | max restarts | n | ended before the budget | never accepted | ended early given never accepted | ended early given accepted | median evaluations |
|---|---|---|---|---|---|---|---|
| 100 | 0 | 900 | 174/900 (19%, 17-22) | 158/900 (18%, 15-20) | 157/158 (99%, 97-100) | 17/742 (2%, 1-4) | 101 |
| 100 | 2 | 900 | 129/900 (14%, 12-17) | 101/900 (11%, 9-13) | 101/101 (100%, 96-100) | 28/799 (4%, 2-5) | 101 |
| 100 | 5 | 900 | 96/900 (11%, 9-13) | 61/900 (7%, 5-9) | 59/61 (97%, 89-99) | 37/839 (4%, 3-6) | 101 |
| 500 | 0 | 900 | 180/900 (20%, 18-23) | 157/900 (17%, 15-20) | 157/157 (100%, 98-100) | 23/743 (3%, 2-5) | 501 |
| 500 | 2 | 900 | 139/900 (15%, 13-18) | 106/900 (12%, 10-14) | 106/106 (100%, 97-100) | 33/794 (4%, 3-6) | 501 |
| 500 | 5 | 900 | 106/900 (12%, 10-14) | 73/900 (8%, 7-10) | 73/73 (100%, 95-100) | 33/827 (4%, 3-6) | 501 |
| 1000 | 0 | 900 | 179/900 (20%, 17-23) | 157/900 (17%, 15-20) | 157/157 (100%, 98-100) | 22/743 (3%, 2-4) | 1001 |
| 1000 | 2 | 900 | 135/900 (15%, 13-17) | 101/900 (11%, 9-13) | 101/101 (100%, 96-100) | 34/799 (4%, 3-6) | 1001 |
| 1000 | 5 | 900 | 108/900 (12%, 10-14) | 68/900 (8%, 6-9) | 68/68 (100%, 95-100) | 40/832 (5%, 4-6) | 1001 |
| 3500 | 0 | 900 | 175/900 (19%, 17-22) | 150/900 (17%, 14-19) | 150/150 (100%, 98-100) | 25/750 (3%, 2-5) | 3501 |
| 3500 | 2 | 900 | 134/900 (15%, 13-17) | 101/900 (11%, 9-13) | 101/101 (100%, 96-100) | 33/799 (4%, 3-6) | 3501 |
| 3500 | 5 | 900 | 106/900 (12%, 10-14) | 68/900 (8%, 6-9) | 68/68 (100%, 95-100) | 38/832 (5%, 3-6) | 3501 |
<!-- /table:mc_stopping -->

<!-- table:mg_budget -->
| optimizer | budget class | config | evaluations | < 1 deg, 1 deg start | < 1 deg, 2 deg start | < 0.5 deg, 1 deg start | median error of the < 1 deg runs | < 0.02 deg, 1 deg start | all configs in class, < 1 deg | all configs in class, < 0.5 deg |
|---|---|---|---|---|---|---|---|---|---|---|
| riemannian_adam_geoopt | <=200 steps | lr 0.0005, 100 steps, b1 0.95 | 101 | 88/100 (88%, 80-93) | 58/100 (58%, 48-67) | 61/100 (61%, 51-70) | 0.330 | 2/100 (2%, 1-7) | 1339/2800 (48%, 46-50) | 730/2800 (26%, 24-28) |
| riemannian_adam_geoopt | 500-1000 steps | lr 0.0001, 500 steps, b1 0.9 | 501 | 88/100 (88%, 80-93) | 57/100 (57%, 47-66) | 61/100 (61%, 51-70) | 0.352 | 1/100 (1%, 0-5) | 533/1400 (38%, 36-41) | 274/1400 (20%, 18-22) |
| riemannian_adam_manual | <=200 steps | lr 0.0001, 200 steps, b1 0.9 | 201 | 94/100 (94%, 88-97) | 49/100 (49%, 39-59) | 58/100 (58%, 48-67) | 0.415 | 5/100 (5%, 2-11) | 1287/2800 (46%, 44-48) | 680/2800 (24%, 23-26) |
| riemannian_adam_manual | 500-1000 steps | lr 0.0001, 500 steps, b1 0.9 | 501 | 86/100 (86%, 78-91) | 53/100 (53%, 43-62) | 53/100 (53%, 43-62) | 0.432 | 9/100 (9%, 5-16) | 514/1400 (37%, 34-39) | 271/1400 (19%, 17-22) |
| riemannian_sgd_plain | <=200 steps | lr 0.0001, 100 steps | 101 | 88/100 (88%, 80-93) | 52/100 (52%, 42-62) | 54/100 (54%, 44-63) | 0.420 | 4/100 (4%, 2-10) | 292/1000 (29%, 26-32) | 135/1000 (14%, 12-16) |
| riemannian_sgd_plain | 500-1000 steps | lr 0.0001, 500 steps | 501 | 56/100 (56%, 46-65) | 37/100 (37%, 28-47) | 26/100 (26%, 18-35) | 0.534 | 3/100 (3%, 1-8) | 110/500 (22%, 19-26) | 46/500 (9%, 7-12) |
| riemannian_sgd_momentum | <=200 steps | lr 1e-05, 200 steps, mom 0.5 | 201 | 94/100 (94%, 88-97) | 51/100 (51%, 41-61) | 59/100 (59%, 49-68) | 0.357 | 1/100 (1%, 0-5) | 567/1500 (38%, 35-40) | 272/1500 (18%, 16-20) |
| riemannian_sgld | <=200 steps | lr 0.0001, 100 steps, T 0.001 | 101 | 92/100 (92%, 85-96) | 50/100 (50%, 40-60) | 62/100 (62%, 52-71) | 0.402 | 0/100 (0%, 0-4) | 842/3000 (28%, 26-30) | 459/3000 (15%, 14-17) |
| riemannian_sgld | 500-1000 steps | lr 0.0001, 500 steps, T 0.001 | 501 | 75/100 (75%, 66-82) | 58/100 (58%, 48-67) | 47/100 (47%, 38-57) | 0.441 | 0/100 (0%, 0-4) | 273/1500 (18%, 16-20) | 143/1500 (10%, 8-11) |
| mc_optimizer | <=200 steps | 100 steps, 2 restarts, step frac 0.5 | 101 | 88/100 (88%, 80-93) | 36/100 (36%, 27-46) | 49/100 (49%, 39-59) | 0.465 | 2/100 (2%, 1-7) | 635/900 (71%, 67-73) | 242/900 (27%, 24-30) |
| mc_optimizer | 500-1000 steps | 500 steps, 2 restarts, step frac 0.5 | 501 | 88/100 (88%, 80-93) | 32/100 (32%, 24-42) | 49/100 (49%, 39-59) | 0.464 | 1/100 (1%, 0-5) | 1299/1800 (72%, 70-74) | 479/1800 (27%, 25-29) |
| mc_optimizer | 3500 steps | 3500 steps, 5 restarts, step frac 0.5 | 3501 | 90/100 (90%, 83-94) | 35/100 (35%, 26-45) | 50/100 (50%, 40-60) | 0.461 | 2/100 (2%, 1-7) | 668/900 (74%, 71-77) | 245/900 (27%, 24-30) |
<!-- /table:mg_budget -->

<!-- table:mg_by_rperp -->
| runs | r_perp | < 1 deg | < 0.5 deg | median error of the < 1 deg runs |
|---|---|---|---|---|
| Adam (geoopt), recorded config | <277 um | 31/33 (94%, 80-98) | 17/33 (52%, 35-67) | 0.468 |
| Adam (geoopt), recorded config | 277-412 um | 33/33 (100%, 90-100) | 18/33 (55%, 38-70) | 0.421 |
| Adam (geoopt), recorded config | >=412 um | 32/34 (94%, 81-98) | 18/34 (53%, 37-69) | 0.483 |
| MC, recorded config | <277 um | 29/33 (88%, 73-95) | 14/33 (42%, 27-59) | 0.517 |
| MC, recorded config | 277-412 um | 29/33 (88%, 73-95) | 13/33 (39%, 25-56) | 0.531 |
| MC, recorded config | >=412 um | 33/34 (97%, 85-99) | 12/34 (35%, 21-52) | 0.576 |
| Adam (geoopt), all configs | <277 um | 624/1386 (45%, 42-48) | 305/1386 (22%, 20-24) | 0.511 |
| Adam (geoopt), all configs | 277-412 um | 679/1386 (49%, 46-52) | 370/1386 (27%, 24-29) | 0.446 |
| Adam (geoopt), all configs | >=412 um | 569/1428 (40%, 37-42) | 329/1428 (23%, 21-25) | 0.426 |
| MC, all configs | <277 um | 865/1188 (73%, 70-75) | 318/1188 (27%, 24-29) | 0.594 |
| MC, all configs | 277-412 um | 861/1188 (72%, 70-75) | 326/1188 (27%, 25-30) | 0.590 |
| MC, all configs | >=412 um | 876/1224 (72%, 69-74) | 322/1224 (26%, 24-29) | 0.596 |
| Adam (geoopt), recorded config (same hp id) | 12 um (ThreeVoxels) | 3/3 | 1/3 | 0.569 |
| MC, recorded config (same hp id) | 12 um (ThreeVoxels) | 3/3 | 1/3 | 0.791 |
| Adam (geoopt), all configs (3 voxels x configs) | 12 um (ThreeVoxels) | 79/126 (63%, 54-71) | 38/126 (30%, 23-39) | 0.526 |
| MC, all configs (3 voxels x configs) | 12 um (ThreeVoxels) | 62/108 (57%, 48-66) | 14/108 (13%, 8-21) | 0.732 |
<!-- /table:mg_by_rperp -->

<!-- table:mg_headline -->
| optimizer | config | evaluations | < 1 deg, 1 deg start | < 1 deg, 2 deg start | < 1 deg, 5 deg start | < 0.5 deg, 1 deg start | < 0.5 deg, 2 deg start | recorded % |
|---|---|---|---|---|---|---|---|---|
| riemannian_adam_geoopt | lr 0.0001, 100 steps, b1 0.9 | 101 | 96/100 (96%, 90-98) | 41/100 (41%, 32-51) | 0/100 (0%, 0-4) | 53/100 (53%, 43-62) | 0/100 (0%, 0-4) | 96/41/0 |
| riemannian_adam_manual | lr 0.0001, 200 steps, b1 0.9 | 201 | 94/100 (94%, 88-97) | 49/100 (49%, 39-59) | 0/100 (0%, 0-4) | 58/100 (58%, 48-67) | 27/100 (27%, 19-36) | 94/49/0 |
| riemannian_sgd_plain | lr 0.0001, 100 steps | 101 | 88/100 (88%, 80-93) | 52/100 (52%, 42-62) | 0/100 (0%, 0-4) | 54/100 (54%, 44-63) | 25/100 (25%, 18-34) | 88/52/0 |
| riemannian_sgd_momentum | lr 1e-05, 200 steps, mom 0.5 | 201 | 94/100 (94%, 88-97) | 51/100 (51%, 41-61) | 0/100 (0%, 0-4) | 59/100 (59%, 49-68) | 21/100 (21%, 14-30) | 94/51/0 |
| riemannian_sgld | lr 0.0001, 100 steps, T 0.001 | 101 | 92/100 (92%, 85-96) | 50/100 (50%, 40-60) | 0/100 (0%, 0-4) | 62/100 (62%, 52-71) | 18/100 (18%, 12-27) | 92/50/0 |
| mc_optimizer | 3500 steps, 2 restarts, step frac 0.5 | 3501 | 91/100 (91%, 84-95) | 40/100 (40%, 31-50) | 6/100 (6%, 3-12) | 39/100 (39%, 30-49) | 18/100 (18%, 12-27) | 92/40/6 |
<!-- /table:mg_headline -->

<!-- table:mg_precision -->
| optimizer | runs | n | n < 1 deg | error of the < 1 deg runs q25/50/75/90 | < 0.5 deg | < 0.1 deg | < 0.05 deg | < 0.02 deg |
|---|---|---|---|---|---|---|---|---|
| riemannian_adam_geoopt | recorded config | 100 | 96 | 0.262/0.472/0.705/0.857 | 53/100 (53%, 43-62) | 3/100 (3%, 1-8) | 1/100 (1%, 0-5) | 0/100 (0%, 0-4) |
| riemannian_adam_geoopt | all configs | 4200 | 1872 | 0.261/0.468/0.717/0.867 | 1004/4200 (24%, 23-25) | 129/4200 (3%, 3-4) | 57/4200 (1%, 1-2) | 47/4200 (1%, 1-1) |
| riemannian_adam_manual | recorded config | 100 | 94 | 0.245/0.415/0.647/0.840 | 58/100 (58%, 48-67) | 12/100 (12%, 7-20) | 5/100 (5%, 2-11) | 5/100 (5%, 2-11) |
| riemannian_adam_manual | all configs | 4200 | 1801 | 0.270/0.473/0.734/0.890 | 951/4200 (23%, 21-24) | 141/4200 (3%, 3-4) | 79/4200 (2%, 2-2) | 62/4200 (1%, 1-2) |
| riemannian_sgd_plain | recorded config | 100 | 88 | 0.229/0.420/0.657/0.854 | 54/100 (54%, 44-63) | 10/100 (10%, 6-17) | 5/100 (5%, 2-11) | 4/100 (4%, 2-10) |
| riemannian_sgd_plain | all configs | 1500 | 402 | 0.347/0.542/0.795/0.920 | 181/1500 (12%, 11-14) | 23/1500 (2%, 1-2) | 14/1500 (1%, 1-2) | 12/1500 (1%, 0-1) |
| riemannian_sgd_momentum | recorded config | 100 | 94 | 0.212/0.357/0.658/0.877 | 59/100 (59%, 49-68) | 11/100 (11%, 6-19) | 3/100 (3%, 1-8) | 1/100 (1%, 0-5) |
| riemannian_sgd_momentum | all configs | 1500 | 567 | 0.296/0.511/0.749/0.921 | 272/1500 (18%, 16-20) | 38/1500 (3%, 2-3) | 17/1500 (1%, 1-2) | 6/1500 (0%, 0-1) |
| riemannian_sgld | recorded config | 100 | 92 | 0.206/0.402/0.654/0.852 | 62/100 (62%, 52-71) | 4/100 (4%, 2-10) | 1/100 (1%, 0-5) | 0/100 (0%, 0-4) |
| riemannian_sgld | all configs | 4500 | 1115 | 0.250/0.463/0.707/0.870 | 602/4500 (13%, 12-14) | 59/4500 (1%, 1-2) | 6/4500 (0%, 0-0) | 1/4500 (0%, 0-0) |
| mc_optimizer | recorded config | 100 | 91 | 0.312/0.531/0.699/0.839 | 39/100 (39%, 30-49) | 5/100 (5%, 2-11) | 1/100 (1%, 0-5) | 0/100 (0%, 0-4) |
| mc_optimizer | all configs | 3600 | 2602 | 0.378/0.592/0.763/0.882 | 966/3600 (27%, 25-28) | 125/3600 (3%, 3-4) | 61/3600 (2%, 1-2) | 23/3600 (1%, 0-1) |
<!-- /table:mg_precision -->
