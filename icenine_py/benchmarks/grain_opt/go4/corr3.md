## go4_corr3_As100

args: `{"name": "go4_corr3_As100", "pipeline": "Astrict/100", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 1.1181 | 6.2648 | 1020 of 1043 (96.7 to 98.5%) | 553 of 1043 (50.0 to 56.0%) | 0.6292 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.6795 | 5.1923 | 83 of 86 (90.2 to 98.8%) | 38 of 86 (34.2 to 54.7%) | 0.5299 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 1.1632 | 6.3059 | 937 of 957 (96.8 to 98.6%) | 515 of 957 (50.6 to 57.0%) | 0.6381 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.3406 | 0.6701 | 254 of 269 (91.0 to 96.6%) | 2 of 269 (0.2 to 2.7%) | 0.3080 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.6523 | 1.7809 | 271 of 277 (95.4 to 99.0%) | 93 of 277 (28.3 to 39.3%) | 0.5269 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 1.6857 | 3.3761 | 243 of 244 (97.7 to 99.9%) | 208 of 244 (80.3 to 89.1%) | 0.7995 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 4.4092 | 8.2467 | 252 of 253 (97.8 to 99.9%) | 250 of 253 (96.6 to 99.6%) | 0.9186 |

grains: 20; stages Astrict/100; accepted per stage [20]; objective evaluations per voxel 113.0; stage wall time 3.0 s (2.9 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.37 > 1.5); single run, one worker

## go4_corr3_As100_V61

args: `{"name": "go4_corr3_As100_V61", "pipeline": "Astrict/100,V/61", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 1.1181 | 6.2648 | 1020 of 1043 (96.7 to 98.5%) | 553 of 1043 (50.0 to 56.0%) | 0.6292 |
| all | s1 | 1043 | 1.0912 | 6.2505 | 932 of 1043 (87.3 to 91.1%) | 540 of 1043 (48.7 to 54.8%) | 0.5619 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.6795 | 5.1923 | 83 of 86 (90.2 to 98.8%) | 38 of 86 (34.2 to 54.7%) | 0.5299 |
| interior (all) | s1 | 86 | 0.5727 | 5.1597 | 75 of 86 (78.5 to 92.7%) | 37 of 86 (33.1 to 53.6%) | 0.4607 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 1.1632 | 6.3059 | 937 of 957 (96.8 to 98.6%) | 515 of 957 (50.6 to 57.0%) | 0.6381 |
| boundary (all) | s1 | 957 | 1.1331 | 6.2887 | 857 of 957 (87.5 to 91.3%) | 503 of 957 (49.4 to 55.7%) | 0.5710 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.3406 | 0.6701 | 254 of 269 (91.0 to 96.6%) | 2 of 269 (0.2 to 2.7%) | 0.3080 |
| level 0.5 (all) | s1 | 269 | 0.2236 | 0.6284 | 186 of 269 (63.4 to 74.4%) | 2 of 269 (0.2 to 2.7%) | 0.2018 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.6523 | 1.7809 | 271 of 277 (95.4 to 99.0%) | 93 of 277 (28.3 to 39.3%) | 0.5269 |
| level 1.0 (all) | s1 | 277 | 0.5535 | 1.7592 | 253 of 277 (87.4 to 94.1%) | 86 of 277 (25.9 to 36.7%) | 0.4291 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 1.6857 | 3.3761 | 243 of 244 (97.7 to 99.9%) | 208 of 244 (80.3 to 89.1%) | 0.7995 |
| level 2.0 (all) | s1 | 244 | 1.7058 | 3.4140 | 241 of 244 (96.4 to 99.6%) | 201 of 244 (77.1 to 86.6%) | 0.7609 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 4.4092 | 8.2467 | 252 of 253 (97.8 to 99.9%) | 250 of 253 (96.6 to 99.6%) | 0.9186 |
| level 5.0 (all) | s1 | 253 | 4.3744 | 8.2645 | 252 of 253 (97.8 to 99.9%) | 251 of 253 (97.2 to 99.8%) | 0.8982 |

grains: 20; stages Astrict/100,V/61; accepted per stage [20, 20]; objective evaluations per voxel 177.0; stage wall time 8.9 s (8.5 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.58 > 1.5); single run, one worker

## go4_corr3_As100_V61s03

args: `{"name": "go4_corr3_As100_V61s03", "pipeline": "Astrict/100,V/61/0.3", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 1.1181 | 6.2648 | 1020 of 1043 (96.7 to 98.5%) | 553 of 1043 (50.0 to 56.0%) | 0.6292 |
| all | s1 | 1043 | 0.3732 | 6.3375 | 731 of 1043 (67.2 to 72.8%) | 387 of 1043 (34.2 to 40.1%) | 0.4432 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.6795 | 5.1923 | 83 of 86 (90.2 to 98.8%) | 38 of 86 (34.2 to 54.7%) | 0.5299 |
| interior (all) | s1 | 86 | 0.1878 | 4.3786 | 51 of 86 (48.7 to 69.1%) | 23 of 86 (18.5 to 36.9%) | 0.3360 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 1.1632 | 6.3059 | 937 of 957 (96.8 to 98.6%) | 515 of 957 (50.6 to 57.0%) | 0.6381 |
| boundary (all) | s1 | 957 | 0.4003 | 6.4323 | 680 of 957 (68.1 to 73.8%) | 364 of 957 (35.0 to 41.2%) | 0.4528 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.3406 | 0.6701 | 254 of 269 (91.0 to 96.6%) | 2 of 269 (0.2 to 2.7%) | 0.3080 |
| level 0.5 (all) | s1 | 269 | 0.0748 | 0.2856 | 87 of 269 (27.0 to 38.1%) | 0 of 269 (0.0 to 1.4%) | 0.1345 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.6523 | 1.7809 | 271 of 277 (95.4 to 99.0%) | 93 of 277 (28.3 to 39.3%) | 0.5269 |
| level 1.0 (all) | s1 | 277 | 0.1538 | 0.9646 | 164 of 277 (53.3 to 64.8%) | 13 of 277 (2.8 to 7.9%) | 0.2486 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 1.6857 | 3.3761 | 243 of 244 (97.7 to 99.9%) | 208 of 244 (80.3 to 89.1%) | 0.7995 |
| level 2.0 (all) | s1 | 244 | 1.0545 | 3.5481 | 228 of 244 (89.6 to 95.9%) | 130 of 244 (47.0 to 59.4%) | 0.5804 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 4.4092 | 8.2467 | 252 of 253 (97.8 to 99.9%) | 250 of 253 (96.6 to 99.6%) | 0.9186 |
| level 5.0 (all) | s1 | 253 | 4.2614 | 8.4879 | 252 of 253 (97.8 to 99.9%) | 244 of 253 (93.4 to 98.1%) | 0.8522 |

grains: 20; stages Astrict/100,V/61/0.3; accepted per stage [20, 20]; objective evaluations per voxel 177.0; stage wall time 9.5 s (9.1 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.11 > 1.5); single run, one worker

## go4_corr3_V61s03

args: `{"name": "go4_corr3_V61s03", "pipeline": "V/61/0.3", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 0.3648 | 6.0194 | 777 of 1043 (71.8 to 77.0%) | 378 of 1043 (33.4 to 39.2%) | 0.4485 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.1507 | 5.2364 | 56 of 86 (54.6 to 74.3%) | 20 of 86 (15.6 to 33.2%) | 0.3283 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 0.4063 | 6.2075 | 721 of 957 (72.5 to 78.0%) | 358 of 957 (34.4 to 40.5%) | 0.4593 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.0796 | 0.2991 | 96 of 269 (30.2 to 41.6%) | 0 of 269 (0.0 to 1.4%) | 0.1395 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.1826 | 1.0182 | 197 of 277 (65.5 to 76.1%) | 15 of 277 (3.3 to 8.7%) | 0.2589 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 0.9593 | 3.1135 | 232 of 244 (91.6 to 97.2%) | 118 of 244 (42.2 to 54.6%) | 0.5772 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 4.2074 | 8.0517 | 252 of 253 (97.8 to 99.9%) | 245 of 253 (93.9 to 98.4%) | 0.8605 |

grains: 20; stages V/61/0.3; accepted per stage [20]; objective evaluations per voxel 64.0; stage wall time 6.7 s (6.4 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.46 > 1.5); single run, one worker

## go4_corr3_Als150_As100_V61s03

args: `{"name": "go4_corr3_Als150_As100_V61s03", "pipeline": "Als/150,Astrict/100,V/61/0.3", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 0.9783 | 6.0999 | 1016 of 1043 (96.3 to 98.2%) | 514 of 1043 (46.3 to 52.3%) | 0.6254 |
| all | s1 | 1043 | 0.9129 | 6.0225 | 1010 of 1043 (95.6 to 97.7%) | 497 of 1043 (44.6 to 50.7%) | 0.5927 |
| all | s2 | 1043 | 0.2797 | 6.2922 | 711 of 1043 (65.3 to 70.9%) | 366 of 1043 (32.3 to 38.0%) | 0.4254 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.7963 | 4.8811 | 84 of 86 (91.9 to 99.4%) | 35 of 86 (30.9 to 51.3%) | 0.5497 |
| interior (all) | s1 | 86 | 0.6561 | 5.0596 | 85 of 86 (93.7 to 99.8%) | 30 of 86 (25.7 to 45.4%) | 0.4984 |
| interior (all) | s2 | 86 | 0.1404 | 4.1229 | 51 of 86 (48.7 to 69.1%) | 23 of 86 (18.5 to 36.9%) | 0.3266 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 1.0011 | 6.1403 | 932 of 957 (96.2 to 98.2%) | 479 of 957 (46.9 to 53.2%) | 0.6322 |
| boundary (all) | s1 | 957 | 0.9409 | 6.1721 | 925 of 957 (95.3 to 97.6%) | 467 of 957 (45.6 to 52.0%) | 0.6012 |
| boundary (all) | s2 | 957 | 0.2898 | 6.3405 | 660 of 957 (66.0 to 71.8%) | 343 of 957 (32.9 to 38.9%) | 0.4343 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.2450 | 0.8378 | 248 of 269 (88.4 to 94.8%) | 6 of 269 (1.0 to 4.8%) | 0.3109 |
| level 0.5 (all) | s1 | 269 | 0.2433 | 0.7276 | 249 of 269 (88.8 to 95.1%) | 4 of 269 (0.6 to 3.8%) | 0.2666 |
| level 0.5 (all) | s2 | 269 | 0.0713 | 0.2807 | 87 of 269 (27.0 to 38.1%) | 0 of 269 (0.0 to 1.4%) | 0.1343 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.5716 | 1.4426 | 273 of 277 (96.3 to 99.4%) | 57 of 277 (16.2 to 25.7%) | 0.4879 |
| level 1.0 (all) | s1 | 277 | 0.5521 | 1.3313 | 268 of 277 (93.9 to 98.3%) | 48 of 277 (13.3 to 22.2%) | 0.4519 |
| level 1.0 (all) | s2 | 277 | 0.1083 | 0.6351 | 147 of 277 (47.2 to 58.9%) | 3 of 277 (0.4 to 3.1%) | 0.2106 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 1.6674 | 3.2903 | 242 of 244 (97.1 to 99.8%) | 204 of 244 (78.4 to 87.7%) | 0.8082 |
| level 2.0 (all) | s1 | 244 | 1.6875 | 3.3394 | 240 of 244 (95.9 to 99.4%) | 197 of 244 (75.3 to 85.2%) | 0.7762 |
| level 2.0 (all) | s2 | 244 | 0.9291 | 3.0653 | 226 of 244 (88.6 to 95.3%) | 116 of 244 (41.4 to 53.8%) | 0.5556 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 4.3092 | 8.3307 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9338 |
| level 5.0 (all) | s1 | 253 | 4.2506 | 8.4318 | 253 of 253 (98.5 to 100.0%) | 248 of 253 (95.5 to 99.2%) | 0.9166 |
| level 5.0 (all) | s2 | 253 | 4.1913 | 8.4377 | 251 of 253 (97.2 to 99.8%) | 247 of 253 (94.9 to 98.9%) | 0.8447 |

grains: 20; stages Als/150,Astrict/100,V/61/0.3; accepted per stage [19, 20, 20]; objective evaluations per voxel 340.0; stage wall time 59.1 s (56.7 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.36 > 1.5); single run, one worker

