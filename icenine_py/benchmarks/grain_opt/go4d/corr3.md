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

## go4d_corr3_L366x3_G15

args: `{"name": "go4d_corr3_L366x3_G15", "pipeline": "L/366/3,G/15/0.15", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 0.0137 | 5.1579 | 290 of 1043 (25.2 to 30.6%) | 199 of 1043 (16.8 to 21.6%) | 0.2048 |
| all | s1 | 1043 | 0.0110 | 2.1481 | 133 of 1043 (10.9 to 14.9%) | 83 of 1043 (6.5 to 9.8%) | 0.1177 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.0151 | 1.2136 | 16 of 86 (11.8 to 28.1%) | 7 of 86 (4.0 to 15.9%) | 0.1379 |
| interior (all) | s1 | 86 | 0.0138 | 0.1279 | 7 of 86 (4.0 to 15.9%) | 1 of 86 (0.2 to 6.3%) | 0.1034 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 0.0134 | 5.9503 | 274 of 957 (25.9 to 31.6%) | 192 of 957 (17.6 to 22.7%) | 0.2108 |
| boundary (all) | s1 | 957 | 0.0108 | 2.5820 | 126 of 957 (11.2 to 15.5%) | 82 of 957 (7.0 to 10.5%) | 0.1190 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.0118 | 0.0724 | 9 of 269 (1.8 to 6.2%) | 0 of 269 (0.0 to 1.4%) | 0.0815 |
| level 0.5 (all) | s1 | 269 | 0.0118 | 0.0674 | 7 of 269 (1.3 to 5.3%) | 1 of 269 (0.1 to 2.1%) | 0.0797 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.0094 | 0.3293 | 25 of 277 (6.2 to 13.0%) | 0 of 277 (0.0 to 1.4%) | 0.0944 |
| level 1.0 (all) | s1 | 277 | 0.0090 | 0.0682 | 13 of 277 (2.8 to 7.9%) | 1 of 277 (0.1 to 2.0%) | 0.0860 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 0.0110 | 1.5366 | 75 of 244 (25.3 to 36.8%) | 25 of 244 (7.0 to 14.7%) | 0.1869 |
| level 2.0 (all) | s1 | 244 | 0.0094 | 0.8518 | 22 of 244 (6.0 to 13.3%) | 10 of 244 (2.2 to 7.4%) | 0.1195 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 2.1502 | 11.2938 | 181 of 253 (65.7 to 76.7%) | 174 of 253 (62.8 to 74.2%) | 0.4740 |
| level 5.0 (all) | s1 | 253 | 0.0287 | 46.2195 | 91 of 253 (30.3 to 42.1%) | 71 of 253 (22.9 to 33.9%) | 0.1911 |

grains: 20; stages L/366/3,G/15/0.15; accepted per stage [20, 20]; objective evaluations per voxel 1101.7; stage wall time 297.4 s (285.1 ms per voxel); timing: contended: timing preflight refused (Python at 100% CPU, Python at 100% CPU, Python at 100% CPU); single run, one worker

## go4d_corr3_L366x3_G10

args: `{"name": "go4d_corr3_L366x3_G10", "pipeline": "L/366/3,G/10/0.10", "sample": "corr3", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1421 | 6.0746 | 1039 of 1043 (99.0 to 99.9%) | 574 of 1043 (52.0 to 58.0%) | 0.7105 |
| all | s0 | 1043 | 0.0137 | 5.1579 | 290 of 1043 (25.2 to 30.6%) | 199 of 1043 (16.8 to 21.6%) | 0.2048 |
| all | s1 | 1043 | 0.0106 | 2.1481 | 129 of 1043 (10.5 to 14.5%) | 83 of 1043 (6.5 to 9.8%) | 0.1173 |
| interior (all) | start | 86 | 0.9255 | 4.7694 | 86 of 86 (95.7 to 100.0%) | 40 of 86 (36.3 to 57.0%) | 0.5999 |
| interior (all) | s0 | 86 | 0.0151 | 1.2136 | 16 of 86 (11.8 to 28.1%) | 7 of 86 (4.0 to 15.9%) | 0.1379 |
| interior (all) | s1 | 86 | 0.0125 | 0.1279 | 7 of 86 (4.0 to 15.9%) | 1 of 86 (0.2 to 6.3%) | 0.1024 |
| boundary (all) | start | 957 | 1.1895 | 6.2064 | 953 of 957 (98.9 to 99.8%) | 534 of 957 (52.6 to 58.9%) | 0.7205 |
| boundary (all) | s0 | 957 | 0.0134 | 5.9503 | 274 of 957 (25.9 to 31.6%) | 192 of 957 (17.6 to 22.7%) | 0.2108 |
| boundary (all) | s1 | 957 | 0.0105 | 2.5820 | 122 of 957 (10.8 to 15.0%) | 82 of 957 (7.0 to 10.5%) | 0.1186 |
| level 0.5 (all) | start | 269 | 0.4140 | 0.7958 | 265 of 269 (96.2 to 99.4%) | 6 of 269 (1.0 to 4.8%) | 0.4029 |
| level 0.5 (all) | s0 | 269 | 0.0118 | 0.0724 | 9 of 269 (1.8 to 6.2%) | 0 of 269 (0.0 to 1.4%) | 0.0815 |
| level 0.5 (all) | s1 | 269 | 0.0108 | 0.0563 | 6 of 269 (1.0 to 4.8%) | 1 of 269 (0.1 to 2.1%) | 0.0793 |
| level 1.0 (all) | start | 277 | 0.8781 | 1.6497 | 277 of 277 (98.6 to 100.0%) | 106 of 277 (32.7 to 44.1%) | 0.6558 |
| level 1.0 (all) | s0 | 277 | 0.0094 | 0.3293 | 25 of 277 (6.2 to 13.0%) | 0 of 277 (0.0 to 1.4%) | 0.0944 |
| level 1.0 (all) | s1 | 277 | 0.0090 | 0.0565 | 10 of 277 (2.0 to 6.5%) | 1 of 277 (0.1 to 2.0%) | 0.0853 |
| level 2.0 (all) | start | 244 | 1.7731 | 3.1407 | 244 of 244 (98.4 to 100.0%) | 215 of 244 (83.5 to 91.6%) | 0.8651 |
| level 2.0 (all) | s0 | 244 | 0.0110 | 1.5366 | 75 of 244 (25.3 to 36.8%) | 25 of 244 (7.0 to 14.7%) | 0.1869 |
| level 2.0 (all) | s1 | 244 | 0.0091 | 0.8518 | 22 of 244 (6.0 to 13.3%) | 10 of 244 (2.2 to 7.4%) | 0.1190 |
| level 5.0 (all) | start | 253 | 4.3475 | 8.2826 | 253 of 253 (98.5 to 100.0%) | 247 of 253 (94.9 to 98.9%) | 0.9485 |
| level 5.0 (all) | s0 | 253 | 2.1502 | 11.2938 | 181 of 253 (65.7 to 76.7%) | 174 of 253 (62.8 to 74.2%) | 0.4740 |
| level 5.0 (all) | s1 | 253 | 0.0287 | 46.2195 | 91 of 253 (30.3 to 42.1%) | 71 of 253 (22.9 to 33.9%) | 0.1910 |

grains: 20; stages L/366/3,G/10/0.10; accepted per stage [20, 20]; objective evaluations per voxel 1543.0; stage wall time 388.6 s (372.6 ms per voxel); timing: contended: timing preflight refused (Python at 99% CPU, Python at 99% CPU); single run, one worker

