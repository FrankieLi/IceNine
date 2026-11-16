## go4_rand_As100

args: `{"name": "go4_rand_As100", "pipeline": "Astrict/100", "sample": "rand", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 1.2070 | 6.1271 | 1028 of 1043 (97.6 to 99.1%) | 571 of 1043 (51.7 to 57.7%) | 0.6177 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 0.9594 | 5.4135 | 86 of 86 (95.7 to 100.0%) | 42 of 86 (38.6 to 59.2%) | 0.5182 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 1.2212 | 6.1353 | 942 of 957 (97.4 to 99.0%) | 529 of 957 (52.1 to 58.4%) | 0.6266 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.4596 | 0.8748 | 260 of 269 (93.8 to 98.2%) | 2 of 269 (0.2 to 2.7%) | 0.2716 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.8675 | 1.7465 | 274 of 277 (96.9 to 99.6%) | 105 of 277 (32.4 to 43.7%) | 0.5181 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8598 | 3.2018 | 241 of 244 (96.4 to 99.6%) | 213 of 244 (82.5 to 90.9%) | 0.7820 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 4.4833 | 7.7565 | 253 of 253 (98.5 to 100.0%) | 251 of 253 (97.2 to 99.8%) | 0.9362 |

grains: 20; stages Astrict/100; accepted per stage [20]; objective evaluations per voxel 113.0; stage wall time 3.0 s (2.8 ms per voxel); timing: contended: timing preflight refused (Google at 24% CPU, Google at 39% CPU); single run, one worker

## go4_rand_As100_V61

args: `{"name": "go4_rand_As100_V61", "pipeline": "Astrict/100,V/61", "sample": "rand", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 1.2070 | 6.1271 | 1028 of 1043 (97.6 to 99.1%) | 571 of 1043 (51.7 to 57.7%) | 0.6177 |
| all | s1 | 1043 | 1.2117 | 6.1476 | 997 of 1043 (94.2 to 96.7%) | 571 of 1043 (51.7 to 57.7%) | 0.5733 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 0.9594 | 5.4135 | 86 of 86 (95.7 to 100.0%) | 42 of 86 (38.6 to 59.2%) | 0.5182 |
| interior (all) | s1 | 86 | 0.9994 | 5.3841 | 84 of 86 (91.9 to 99.4%) | 43 of 86 (39.7 to 60.3%) | 0.4797 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 1.2212 | 6.1353 | 942 of 957 (97.4 to 99.0%) | 529 of 957 (52.1 to 58.4%) | 0.6266 |
| boundary (all) | s1 | 957 | 1.2193 | 6.1531 | 913 of 957 (93.9 to 96.6%) | 528 of 957 (52.0 to 58.3%) | 0.5817 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.4596 | 0.8748 | 260 of 269 (93.8 to 98.2%) | 2 of 269 (0.2 to 2.7%) | 0.2716 |
| level 0.5 (all) | s1 | 269 | 0.4516 | 0.8466 | 242 of 269 (85.8 to 93.0%) | 3 of 269 (0.4 to 3.2%) | 0.2212 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.8675 | 1.7465 | 274 of 277 (96.9 to 99.6%) | 105 of 277 (32.4 to 43.7%) | 0.5181 |
| level 1.0 (all) | s1 | 277 | 0.8634 | 1.7199 | 263 of 277 (91.7 to 97.0%) | 103 of 277 (31.7 to 43.0%) | 0.4531 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8598 | 3.2018 | 241 of 244 (96.4 to 99.6%) | 213 of 244 (82.5 to 90.9%) | 0.7820 |
| level 2.0 (all) | s1 | 244 | 1.8768 | 3.2847 | 239 of 244 (95.3 to 99.1%) | 213 of 244 (82.5 to 90.9%) | 0.7364 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 4.4833 | 7.7565 | 253 of 253 (98.5 to 100.0%) | 251 of 253 (97.2 to 99.8%) | 0.9362 |
| level 5.0 (all) | s1 | 253 | 4.4746 | 7.7689 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9221 |

grains: 20; stages Astrict/100,V/61; accepted per stage [20, 20]; objective evaluations per voxel 177.0; stage wall time 8.9 s (8.5 ms per voxel); timing: contended: timing preflight refused (Google at 36% CPU); single run, one worker

## go4_rand_As100_V61s03

args: `{"name": "go4_rand_As100_V61s03", "pipeline": "Astrict/100,V/61/0.3", "sample": "rand", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 1.2070 | 6.1271 | 1028 of 1043 (97.6 to 99.1%) | 571 of 1043 (51.7 to 57.7%) | 0.6177 |
| all | s1 | 1043 | 1.1723 | 6.3038 | 910 of 1043 (85.1 to 89.1%) | 562 of 1043 (50.8 to 56.9%) | 0.5100 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 0.9594 | 5.4135 | 86 of 86 (95.7 to 100.0%) | 42 of 86 (38.6 to 59.2%) | 0.5182 |
| interior (all) | s1 | 86 | 0.9716 | 5.5016 | 75 of 86 (78.5 to 92.7%) | 43 of 86 (39.7 to 60.3%) | 0.4214 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 1.2212 | 6.1353 | 942 of 957 (97.4 to 99.0%) | 529 of 957 (52.1 to 58.4%) | 0.6266 |
| boundary (all) | s1 | 957 | 1.2041 | 6.3206 | 835 of 957 (85.0 to 89.2%) | 519 of 957 (51.1 to 57.4%) | 0.5179 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.4596 | 0.8748 | 260 of 269 (93.8 to 98.2%) | 2 of 269 (0.2 to 2.7%) | 0.2716 |
| level 0.5 (all) | s1 | 269 | 0.3138 | 0.8900 | 190 of 269 (64.9 to 75.8%) | 4 of 269 (0.6 to 3.8%) | 0.1932 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.8675 | 1.7465 | 274 of 277 (96.9 to 99.6%) | 105 of 277 (32.4 to 43.7%) | 0.5181 |
| level 1.0 (all) | s1 | 277 | 0.8541 | 1.7603 | 236 of 277 (80.5 to 88.9%) | 110 of 277 (34.1 to 45.6%) | 0.3660 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8598 | 3.2018 | 241 of 244 (96.4 to 99.6%) | 213 of 244 (82.5 to 90.9%) | 0.7820 |
| level 2.0 (all) | s1 | 244 | 1.8616 | 3.6494 | 233 of 244 (92.1 to 97.5%) | 199 of 244 (76.2 to 85.9%) | 0.6406 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 4.4833 | 7.7565 | 253 of 253 (98.5 to 100.0%) | 251 of 253 (97.2 to 99.8%) | 0.9362 |
| level 5.0 (all) | s1 | 253 | 4.5350 | 7.9497 | 251 of 253 (97.2 to 99.8%) | 249 of 253 (96.0 to 99.4%) | 0.8784 |

grains: 20; stages Astrict/100,V/61/0.3; accepted per stage [20, 20]; objective evaluations per voxel 177.0; stage wall time 9.6 s (9.2 ms per voxel); timing: contended: timing preflight refused (Google at 29% CPU, Google at 39% CPU); single run, one worker

## go4_rand_V61s03

args: `{"name": "go4_rand_V61s03", "pipeline": "V/61/0.3", "sample": "rand", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 1.1356 | 6.3296 | 913 of 1043 (85.4 to 89.4%) | 555 of 1043 (50.2 to 56.2%) | 0.5104 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 1.0895 | 5.2606 | 77 of 86 (81.3 to 94.4%) | 44 of 86 (40.8 to 61.4%) | 0.4221 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 1.1502 | 6.3372 | 836 of 957 (85.1 to 89.3%) | 511 of 957 (50.2 to 56.5%) | 0.5183 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.3635 | 0.8872 | 194 of 269 (66.5 to 77.1%) | 5 of 269 (0.8 to 4.3%) | 0.1980 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.8529 | 1.7408 | 235 of 277 (80.1 to 88.6%) | 110 of 277 (34.1 to 45.6%) | 0.3662 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8121 | 3.4284 | 232 of 244 (91.6 to 97.2%) | 194 of 244 (74.0 to 84.1%) | 0.6421 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 4.4762 | 8.0916 | 252 of 253 (97.8 to 99.9%) | 246 of 253 (94.4 to 98.7%) | 0.8734 |

grains: 20; stages V/61/0.3; accepted per stage [20]; objective evaluations per voxel 64.0; stage wall time 6.6 s (6.3 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.58 > 1.5); single run, one worker

## go4_rand_Als150_As100_V61s03

args: `{"name": "go4_rand_Als150_As100_V61s03", "pipeline": "Als/150,Astrict/100,V/61/0.3", "sample": "rand", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 1.1905 | 6.1893 | 1033 of 1043 (98.2 to 99.5%) | 572 of 1043 (51.8 to 57.8%) | 0.6417 |
| all | s1 | 1043 | 1.1908 | 6.2203 | 1023 of 1043 (97.1 to 98.8%) | 573 of 1043 (51.9 to 57.9%) | 0.6112 |
| all | s2 | 1043 | 1.1811 | 6.3103 | 904 of 1043 (84.5 to 88.6%) | 564 of 1043 (51.0 to 57.1%) | 0.4988 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 0.9981 | 5.4081 | 86 of 86 (95.7 to 100.0%) | 43 of 86 (39.7 to 60.3%) | 0.5326 |
| interior (all) | s1 | 86 | 1.0392 | 5.5723 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5206 |
| interior (all) | s2 | 86 | 0.9969 | 5.1549 | 76 of 86 (79.9 to 93.6%) | 43 of 86 (39.7 to 60.3%) | 0.4177 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 1.2161 | 6.2178 | 947 of 957 (98.1 to 99.4%) | 529 of 957 (52.1 to 58.4%) | 0.6515 |
| boundary (all) | s1 | 957 | 1.2065 | 6.2331 | 937 of 957 (96.8 to 98.6%) | 529 of 957 (52.1 to 58.4%) | 0.6194 |
| boundary (all) | s2 | 957 | 1.1990 | 6.3620 | 828 of 957 (84.2 to 88.5%) | 521 of 957 (51.3 to 57.6%) | 0.5060 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.4368 | 0.8352 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.2987 |
| level 0.5 (all) | s1 | 269 | 0.4541 | 0.8731 | 259 of 269 (93.3 to 98.0%) | 1 of 269 (0.1 to 2.1%) | 0.2730 |
| level 0.5 (all) | s2 | 269 | 0.3340 | 0.8646 | 182 of 269 (61.9 to 73.0%) | 3 of 269 (0.4 to 3.2%) | 0.1926 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.8841 | 1.6871 | 276 of 277 (98.0 to 99.9%) | 111 of 277 (34.5 to 45.9%) | 0.5672 |
| level 1.0 (all) | s1 | 277 | 0.8785 | 1.7259 | 271 of 277 (95.4 to 99.0%) | 109 of 277 (33.8 to 45.2%) | 0.5273 |
| level 1.0 (all) | s2 | 277 | 0.8153 | 1.7184 | 238 of 277 (81.3 to 89.5%) | 111 of 277 (34.5 to 45.9%) | 0.3656 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8425 | 3.2849 | 243 of 244 (97.7 to 99.9%) | 210 of 244 (81.2 to 89.9%) | 0.7959 |
| level 2.0 (all) | s1 | 244 | 1.8572 | 3.5445 | 241 of 244 (96.4 to 99.6%) | 213 of 244 (82.5 to 90.9%) | 0.7562 |
| level 2.0 (all) | s2 | 244 | 1.9067 | 3.7682 | 236 of 244 (93.7 to 98.3%) | 203 of 244 (78.0 to 87.4%) | 0.6167 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 4.5089 | 7.7477 | 253 of 253 (98.5 to 100.0%) | 251 of 253 (97.2 to 99.8%) | 0.9391 |
| level 5.0 (all) | s1 | 253 | 4.5287 | 7.9818 | 252 of 253 (97.8 to 99.9%) | 250 of 253 (96.6 to 99.6%) | 0.9229 |
| level 5.0 (all) | s2 | 253 | 4.6843 | 8.3108 | 248 of 253 (95.5 to 99.2%) | 247 of 253 (94.9 to 98.9%) | 0.8563 |

grains: 20; stages Als/150,Astrict/100,V/61/0.3; accepted per stage [15, 20, 20]; objective evaluations per voxel 340.0; stage wall time 48.9 s (46.9 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 2.28 > 1.5); single run, one worker

