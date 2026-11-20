## go4c_corr3_L366x3

args: `{"name": "go4c_corr3_L366x3", "pipeline": "L/366/3", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 0.0137 | 5.1579 | 290 of 1043 (25.2 to 30.6%) | 199 of 1043 (16.8 to 21.6%) | 0.2048 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.0151 | 1.2136 | 16 of 86 (11.8 to 28.1%) | 7 of 86 (4.0 to 15.9%) | 0.1379 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 0.0134 | 5.9503 | 274 of 957 (25.9 to 31.6%) | 192 of 957 (17.6 to 22.7%) | 0.2108 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.0118 | 0.0724 | 9 of 269 (1.8 to 6.2%) | 0 of 269 (0.0 to 1.4%) | 0.0815 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.0094 | 0.3293 | 25 of 277 (6.2 to 13.0%) | 0 of 277 (0.0 to 1.4%) | 0.0944 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 0.0110 | 1.5366 | 75 of 244 (25.3 to 36.8%) | 25 of 244 (7.0 to 14.7%) | 0.1869 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 2.1502 | 11.2938 | 181 of 253 (65.7 to 76.7%) | 174 of 253 (62.8 to 74.2%) | 0.4740 |

grains: 20; stages L/366/3; accepted per stage [20]; objective evaluations per voxel 401.0; stage wall time 36.5 s (35.0 ms per voxel); timing: contended: timing preflight refused (Python at 97% CPU); single run, one worker

## go8b_corr3_L366x3_Q3

args: `{"name": "go8b_corr3_L366x3_Q3", "pipeline": "L/366/3,Q/3/0", "sample": "corr3", "stack": null, "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 0.0137 | 5.1579 | 290 of 1043 (25.2 to 30.6%) | 199 of 1043 (16.8 to 21.6%) | 0.2048 |
| all | s1 | 1043 | 0.0137 | 5.1579 | 290 of 1043 (25.2 to 30.6%) | 199 of 1043 (16.8 to 21.6%) | 0.2047 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.0151 | 1.2136 | 16 of 86 (11.8 to 28.1%) | 7 of 86 (4.0 to 15.9%) | 0.1379 |
| interior (all) | s1 | 86 | 0.0147 | 1.2136 | 16 of 86 (11.8 to 28.1%) | 7 of 86 (4.0 to 15.9%) | 0.1377 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 0.0134 | 5.9503 | 274 of 957 (25.9 to 31.6%) | 192 of 957 (17.6 to 22.7%) | 0.2108 |
| boundary (all) | s1 | 957 | 0.0134 | 5.9503 | 274 of 957 (25.9 to 31.6%) | 192 of 957 (17.6 to 22.7%) | 0.2107 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.0118 | 0.0724 | 9 of 269 (1.8 to 6.2%) | 0 of 269 (0.0 to 1.4%) | 0.0815 |
| level 0.5 (all) | s1 | 269 | 0.0117 | 0.0713 | 9 of 269 (1.8 to 6.2%) | 0 of 269 (0.0 to 1.4%) | 0.0814 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.0094 | 0.3293 | 25 of 277 (6.2 to 13.0%) | 0 of 277 (0.0 to 1.4%) | 0.0944 |
| level 1.0 (all) | s1 | 277 | 0.0093 | 0.3293 | 25 of 277 (6.2 to 13.0%) | 0 of 277 (0.0 to 1.4%) | 0.0943 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 0.0110 | 1.5366 | 75 of 244 (25.3 to 36.8%) | 25 of 244 (7.0 to 14.7%) | 0.1869 |
| level 2.0 (all) | s1 | 244 | 0.0108 | 1.5366 | 75 of 244 (25.3 to 36.8%) | 25 of 244 (7.0 to 14.7%) | 0.1868 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 2.1502 | 11.2938 | 181 of 253 (65.7 to 76.7%) | 174 of 253 (62.8 to 74.2%) | 0.4740 |
| level 5.0 (all) | s1 | 253 | 2.1502 | 11.2938 | 181 of 253 (65.7 to 76.7%) | 174 of 253 (62.8 to 74.2%) | 0.4739 |

grains: 20; stages L/366/3,Q/3/0; accepted per stage [20, 20]; objective evaluations per voxel 402.1; stage wall time 79.3 s (76.0 ms per voxel); timing: contended: timing preflight refused (corespotlightd at 198% CPU); single run, one worker

## go8b_corr3_Q3

args: `{"name": "go8b_corr3_Q3", "pipeline": "Q/3/0", "sample": "corr3", "stack": null, "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 1.1390 | 6.0947 | 1039 of 1043 (99.0 to 99.9%) | 568 of 1043 (51.4 to 57.5%) | 0.7030 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5994 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 1.1675 | 6.1860 | 953 of 957 (98.9 to 99.8%) | 528 of 957 (52.0 to 58.3%) | 0.7123 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.4157 | 0.7714 | 265 of 269 (96.2 to 99.4%) | 5 of 269 (0.8 to 4.3%) | 0.3992 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.8624 | 1.6453 | 277 of 277 (98.6 to 100.0%) | 98 of 277 (30.0 to 41.2%) | 0.6439 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 1.7731 | 3.1806 | 244 of 244 (98.4 to 100.0%) | 216 of 244 (83.9 to 91.9%) | 0.8584 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 4.3094 | 8.2674 | 253 of 253 (98.5 to 100.0%) | 249 of 253 (96.0 to 99.4%) | 0.9408 |

grains: 20; stages Q/3/0; accepted per stage [20]; objective evaluations per voxel 1.5; stage wall time 56.0 s (53.7 ms per voxel); timing: contended: timing preflight refused (corespotlightd at 199% CPU, Python at 100% CPU, Python at 100% CPU); single run, one worker

