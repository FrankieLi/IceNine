## go4c_rand_L366x3

args: `{"name": "go4c_rand_L366x3", "pipeline": "L/366/3", "sample": "rand", "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 0.8724 | 7.5064 | 676 of 1043 (61.9 to 67.7%) | 501 of 1043 (45.0 to 51.1%) | 0.3244 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 1.1466 | 6.6120 | 61 of 86 (60.6 to 79.5%) | 44 of 86 (40.8 to 61.4%) | 0.3216 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 0.8623 | 7.5137 | 615 of 957 (61.2 to 67.2%) | 457 of 957 (44.6 to 50.9%) | 0.3247 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.0469 | 0.8838 | 130 of 269 (42.4 to 54.3%) | 8 of 269 (1.5 to 5.8%) | 0.1316 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.5069 | 1.7914 | 148 of 277 (47.5 to 59.2%) | 101 of 277 (31.0 to 42.3%) | 0.2183 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8918 | 3.9088 | 187 of 244 (70.9 to 81.5%) | 182 of 244 (68.8 to 79.6%) | 0.3963 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 5.3482 | 10.5601 | 211 of 253 (78.3 to 87.5%) | 210 of 253 (77.9 to 87.1%) | 0.5764 |

grains: 20; stages L/366/3; accepted per stage [20]; objective evaluations per voxel 401.0; stage wall time 37.2 s (35.6 ms per voxel); timing: contended: timing preflight refused (Python at 100% CPU); single run, one worker

## go8b_rand_L366x3_Q3

args: `{"name": "go8b_rand_L366x3_Q3", "pipeline": "L/366/3,Q/3/0", "sample": "rand", "stack": null, "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 0.8724 | 7.5064 | 676 of 1043 (61.9 to 67.7%) | 501 of 1043 (45.0 to 51.1%) | 0.3244 |
| all | s1 | 1043 | 0.8724 | 7.5064 | 675 of 1043 (61.8 to 67.6%) | 501 of 1043 (45.0 to 51.1%) | 0.3243 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 1.1466 | 6.6120 | 61 of 86 (60.6 to 79.5%) | 44 of 86 (40.8 to 61.4%) | 0.3216 |
| interior (all) | s1 | 86 | 1.1466 | 6.6120 | 61 of 86 (60.6 to 79.5%) | 44 of 86 (40.8 to 61.4%) | 0.3214 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 0.8623 | 7.5137 | 615 of 957 (61.2 to 67.2%) | 457 of 957 (44.6 to 50.9%) | 0.3247 |
| boundary (all) | s1 | 957 | 0.8623 | 7.5137 | 614 of 957 (61.1 to 67.1%) | 457 of 957 (44.6 to 50.9%) | 0.3246 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.0469 | 0.8838 | 130 of 269 (42.4 to 54.3%) | 8 of 269 (1.5 to 5.8%) | 0.1316 |
| level 0.5 (all) | s1 | 269 | 0.0424 | 0.8838 | 129 of 269 (42.1 to 53.9%) | 8 of 269 (1.5 to 5.8%) | 0.1314 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.5069 | 1.7914 | 148 of 277 (47.5 to 59.2%) | 101 of 277 (31.0 to 42.3%) | 0.2183 |
| level 1.0 (all) | s1 | 277 | 0.5069 | 1.7914 | 148 of 277 (47.5 to 59.2%) | 101 of 277 (31.0 to 42.3%) | 0.2182 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8918 | 3.9088 | 187 of 244 (70.9 to 81.5%) | 182 of 244 (68.8 to 79.6%) | 0.3963 |
| level 2.0 (all) | s1 | 244 | 1.8918 | 3.9088 | 187 of 244 (70.9 to 81.5%) | 182 of 244 (68.8 to 79.6%) | 0.3962 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 5.3482 | 10.5601 | 211 of 253 (78.3 to 87.5%) | 210 of 253 (77.9 to 87.1%) | 0.5764 |
| level 5.0 (all) | s1 | 253 | 5.3482 | 10.5593 | 211 of 253 (78.3 to 87.5%) | 210 of 253 (77.9 to 87.1%) | 0.5764 |

grains: 20; stages L/366/3,Q/3/0; accepted per stage [20, 20]; objective evaluations per voxel 402.1; stage wall time 79.8 s (76.5 ms per voxel); timing: contended: timing preflight refused (corespotlightd at 200% CPU); single run, one worker

## go8b_rand_Q3

args: `{"name": "go8b_rand_Q3", "pipeline": "Q/3/0", "sample": "rand", "stack": null, "start": "centre", "grains": 20, "true_grains": null, "group": "true", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1043 | 1.1350 | 6.1741 | 1035 of 1043 (98.5 to 99.6%) | 565 of 1043 (51.1 to 57.2%) | 0.6638 |
| all | s0 | 1043 | 1.1350 | 6.1202 | 1035 of 1043 (98.5 to 99.6%) | 564 of 1043 (51.0 to 57.1%) | 0.6629 |
| interior (all) | start | 86 | 1.0422 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5605 |
| interior (all) | s0 | 86 | 1.0418 | 5.3678 | 86 of 86 (95.7 to 100.0%) | 44 of 86 (40.8 to 61.4%) | 0.5595 |
| boundary (all) | start | 957 | 1.1600 | 6.1828 | 949 of 957 (98.4 to 99.6%) | 521 of 957 (51.3 to 57.6%) | 0.6731 |
| boundary (all) | s0 | 957 | 1.1600 | 6.1344 | 949 of 957 (98.4 to 99.6%) | 520 of 957 (51.2 to 57.5%) | 0.6722 |
| level 0.5 (all) | start | 269 | 0.4435 | 0.8153 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.3006 |
| level 0.5 (all) | s0 | 269 | 0.4435 | 0.8080 | 261 of 269 (94.2 to 98.5%) | 0 of 269 (0.0 to 1.4%) | 0.2999 |
| level 1.0 (all) | start | 277 | 0.8695 | 1.6640 | 277 of 277 (98.6 to 100.0%) | 101 of 277 (31.0 to 42.3%) | 0.5791 |
| level 1.0 (all) | s0 | 277 | 0.8695 | 1.6683 | 277 of 277 (98.6 to 100.0%) | 100 of 277 (30.7 to 41.9%) | 0.5787 |
| level 2.0 (all) | start | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 2.0 (all) | s0 | 244 | 1.8182 | 3.0472 | 244 of 244 (98.4 to 100.0%) | 212 of 244 (82.1 to 90.6%) | 0.8604 |
| level 5.0 (all) | start | 253 | 4.4775 | 7.8565 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9533 |
| level 5.0 (all) | s0 | 253 | 4.4775 | 7.8200 | 253 of 253 (98.5 to 100.0%) | 252 of 253 (97.8 to 99.9%) | 0.9505 |

grains: 20; stages Q/3/0; accepted per stage [20]; objective evaluations per voxel 1.2; stage wall time 46.3 s (44.4 ms per voxel); timing: contended: timing preflight refused (corespotlightd at 107% CPU, Python at 100% CPU, Python at 100% CPU, Python at 100% CPU); single run, one worker

