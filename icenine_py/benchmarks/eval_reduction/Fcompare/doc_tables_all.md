<!-- table:evalFc_e0_cheapest -->
| offset (deg) | share within 1 deg at least | fewest strict-equivalent evaluations (median) | fewest seconds per run (median) |
|---|---|---|---|
| 1 | 100% | F5 expanding shells (83.2) | F5 expanding shells (0.04) |
| 1 | 99% | F5 expanding shells (83.2) | F5 expanding shells (0.04) |
| 1 | 95% | F5 expanding shells (83.2) | F5 expanding shells (0.04) |
| 1 | 90% | F5 expanding shells (83.2) | F5 expanding shells (0.04) |
| 2 | 100% | F5 expanding shells (219) | F5 expanding shells (0.10) |
| 2 | 99% | F5 expanding shells (219) | F5 expanding shells (0.10) |
| 2 | 95% | F5 expanding shells (219) | F5 expanding shells (0.10) |
| 2 | 90% | F5 expanding shells (219) | F5 expanding shells (0.10) |
| 3 | 100% | F4, 5-frame window (331) | Policy, 1-frame window (0.13) |
| 3 | 99% | F4, 5-frame window (331) | Policy, 1-frame window (0.13) |
| 3 | 95% | F4, 5-frame window (331) | Policy, 1-frame window (0.13) |
| 3 | 90% | F4, 5-frame window (331) | Policy, 1-frame window (0.13) |
| 5 | 100% | F4, 5-frame window (370) | Policy, 1-frame window (0.13) |
| 5 | 99% | F4, 5-frame window (370) | Policy, 1-frame window (0.13) |
| 5 | 95% | F4, 1-frame window (358) | Policy, 1-frame window (0.13) |
| 5 | 90% | F4, 1-frame window (358) | Policy, 1-frame window (0.13) |
<!-- /table:evalFc_e0_cheapest -->

<!-- table:evalFc_e0_contention -->
| offset (deg) | option | single-worker, median (s) | contended, same runs, median (s) | ratio of medians |
|---|---|---|---|---|
| 1 | F4, 1-frame window | 0.10 | 0.13 | 1.26 |
| 1 | F4, 5-frame window | 0.11 | 0.14 | 1.22 |
| 1 | F5 expanding shells | 0.04 | 0.04 | 1.15 |
| 1 | F5 then F4 | 0.10 | 0.11 | 1.18 |
| 1 | Policy, 1-frame window | 0.10 | 0.13 | 1.26 |
| 1 | Policy, 5-frame window | 0.11 | 0.14 | 1.22 |
| 1 | CMA-ES, sigma0 1 deg | 0.14 | 0.17 | 1.20 |
| 1 | CMA-ES, sigma0 2 deg | 0.16 | 0.19 | 1.21 |
| 2 | F4, 1-frame window | 0.13 | 0.15 | 1.16 |
| 2 | F4, 5-frame window | 0.13 | 0.15 | 1.16 |
| 2 | F5 expanding shells | 0.10 | 0.11 | 1.13 |
| 2 | F5 then F4 | 0.17 | 0.20 | 1.14 |
| 2 | Policy, 1-frame window | 0.13 | 0.15 | 1.16 |
| 2 | Policy, 5-frame window | 0.13 | 0.15 | 1.16 |
| 2 | CMA-ES, sigma0 1 deg | 0.15 | 0.18 | 1.23 |
| 2 | CMA-ES, sigma0 2 deg | 0.15 | 0.17 | 1.18 |
| 3 | F4, 1-frame window | 0.13 | 0.15 | 1.14 |
| 3 | F4, 5-frame window | 0.14 | 0.16 | 1.18 |
| 3 | F5 expanding shells | 0.22 | 0.24 | 1.11 |
| 3 | F5 then F4 | 0.28 | 0.31 | 1.14 |
| 3 | Policy, 1-frame window | 0.13 | 0.15 | 1.14 |
| 3 | Policy, 5-frame window | 0.14 | 0.16 | 1.18 |
| 3 | CMA-ES, sigma0 1 deg | 0.17 | 0.19 | 1.14 |
| 3 | CMA-ES, sigma0 2 deg | 0.14 | 0.18 | 1.29 |
| 5 | F4, 1-frame window | 0.13 | 0.17 | 1.23 |
| 5 | F4, 5-frame window | 0.15 | 0.18 | 1.19 |
| 5 | F5 expanding shells | 0.58 | 0.68 | 1.16 |
| 5 | F5 then F4 | 0.65 | 0.75 | 1.15 |
| 5 | Policy, 1-frame window | 0.13 | 0.17 | 1.23 |
| 5 | Policy, 5-frame window | 0.15 | 0.18 | 1.19 |
| 5 | CMA-ES, sigma0 1 deg | 0.23 | 0.29 | 1.26 |
| 5 | CMA-ES, sigma0 2 deg | 0.21 | 0.26 | 1.26 |
<!-- /table:evalFc_e0_contention -->

<!-- table:evalFc_e0_main_1deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 377 of 480, 78.5% (74.6-82.0) | n/a | 268 | 342 | 0.10 | 0.15 |
| F4, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 394 of 480, 82.1% (78.4-85.3) | n/a | 274 | 391 | 0.11 | 0.20 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 83.2 | 83.2 | 0.04 | 0.04 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 368 of 480, 76.7% (72.7-80.2) | n/a | 243 | 320 | 0.10 | 0.11 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 377 of 480, 78.5% (74.6-82.0) | 0 of 480 | 268 | 342 | 0.10 | 0.15 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 394 of 480, 82.1% (78.4-85.3) | 0 of 480 | 274 | 391 | 0.11 | 0.20 |
| CMA-ES, sigma0 1 deg | 477 of 480, 99.4% (98.2-99.8) | 476 of 480, 99.2% (97.9-99.7) | n/a | 281 | 358 | 0.14 | 0.18 |
| CMA-ES, sigma0 2 deg | 451 of 480, 94.0% (91.5-95.8) | 442 of 480, 92.1% (89.3-94.2) | n/a | 337 | 526 | 0.16 | 0.35 |
<!-- /table:evalFc_e0_main_1deg -->

<!-- table:evalFc_e0_main_2deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 379 of 480, 79.0% (75.1-82.4) | n/a | 304 | 385 | 0.13 | 0.15 |
| F4, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 392 of 480, 81.7% (78.0-84.9) | n/a | 304 | 429 | 0.13 | 0.21 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 219 | 279 | 0.10 | 0.12 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 372 of 480, 77.5% (73.6-81.0) | n/a | 390 | 474 | 0.17 | 0.32 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 379 of 480, 79.0% (75.1-82.4) | 0 of 480 | 304 | 385 | 0.13 | 0.15 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 392 of 480, 81.7% (78.0-84.9) | 0 of 480 | 304 | 429 | 0.13 | 0.21 |
| CMA-ES, sigma0 1 deg | 450 of 480, 93.8% (91.2-95.6) | 442 of 480, 92.1% (89.3-94.2) | n/a | 323 | 435 | 0.15 | 0.16 |
| CMA-ES, sigma0 2 deg | 430 of 480, 89.6% (86.5-92.0) | 416 of 480, 86.7% (83.3-89.4) | n/a | 344 | 709 | 0.15 | 0.17 |
<!-- /table:evalFc_e0_main_2deg -->

<!-- table:evalFc_e0_main_3deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 476 of 480, 99.2% (97.9-99.7) | 380 of 480, 79.2% (75.3-82.6) | n/a | 331 | 404 | 0.13 | 0.14 |
| F4, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 396 of 480, 82.5% (78.8-85.6) | n/a | 331 | 447 | 0.14 | 0.21 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 35 of 480, 7.3% (5.3-10.0) | n/a | 516 | 516 | 0.22 | 0.23 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 377 of 480, 78.5% (74.6-82.0) | n/a | 648 | 730 | 0.28 | 0.31 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 384 of 480, 80.0% (76.2-83.3) | 4 of 480 | 331 | 412 | 0.13 | 0.14 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 396 of 480, 82.5% (78.8-85.6) | 0 of 480 | 331 | 447 | 0.14 | 0.21 |
| CMA-ES, sigma0 1 deg | 364 of 480, 75.8% (71.8-79.4) | 358 of 480, 74.6% (70.5-78.3) | n/a | 376 | 911 | 0.17 | 0.21 |
| CMA-ES, sigma0 2 deg | 390 of 480, 81.2% (77.5-84.5) | 382 of 480, 79.6% (75.7-82.9) | n/a | 379 | 1109 | 0.14 | 0.87 |
<!-- /table:evalFc_e0_main_3deg -->

<!-- table:evalFc_e0_main_5deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 457 of 480, 95.2% (92.9-96.8) | 369 of 480, 76.9% (72.9-80.4) | n/a | 358 | 454 | 0.13 | 0.15 |
| F4, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 401 of 480, 83.5% (80.0-86.6) | n/a | 370 | 475 | 0.15 | 0.23 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 1364 | 1568 | 0.58 | 0.66 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 373 of 480, 77.7% (73.8-81.2) | n/a | 1533 | 1709 | 0.65 | 0.73 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 388 of 480, 80.8% (77.1-84.1) | 23 of 480 | 377 | 489 | 0.13 | 0.15 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 401 of 480, 83.5% (80.0-86.6) | 0 of 480 | 370 | 475 | 0.15 | 0.23 |
| CMA-ES, sigma0 1 deg | 172 of 480, 35.8% (31.7-40.2) | 165 of 480, 34.4% (30.3-38.7) | n/a | 606 | 1591 | 0.23 | 0.48 |
| CMA-ES, sigma0 2 deg | 261 of 480, 54.4% (49.9-58.8) | 248 of 480, 51.7% (47.2-56.1) | n/a | 554 | 2220 | 0.21 | 0.95 |
<!-- /table:evalFc_e0_main_5deg -->

<!-- table:evalFc_e0_paired -->
| offset (deg) | comparison | fixed | broken | both wrong | both right | voxels with a fixed run | voxels with a broken run |
|---|---|---|---|---|---|---|---|
| 1 | Policy 1 against F4 (1-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 1 | Policy 5 against F4 (5-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 1 | F4 5-frame against F4 1-frame | 0 | 0 | 0 | 480 | 0 | 0 |
| 1 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
| 2 | Policy 1 against F4 (1-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 2 | Policy 5 against F4 (5-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 2 | F4 5-frame against F4 1-frame | 0 | 0 | 0 | 480 | 0 | 0 |
| 2 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
| 3 | Policy 1 against F4 (1-frame) | 4 | 0 | 0 | 476 | 4 | 0 |
| 3 | Policy 5 against F4 (5-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 3 | F4 5-frame against F4 1-frame | 4 | 0 | 0 | 476 | 4 | 0 |
| 3 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
| 5 | Policy 1 against F4 (1-frame) | 23 | 0 | 0 | 457 | 16 | 0 |
| 5 | Policy 5 against F4 (5-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 5 | F4 5-frame against F4 1-frame | 23 | 0 | 0 | 457 | 16 | 0 |
| 5 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
<!-- /table:evalFc_e0_paired -->

<!-- table:evalFc_e0_pervoxel_1deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 121 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 1053 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 1078 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 19 of 20 |
| 1260 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 6674 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 |
| 8214 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 8814 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 11086 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 18 of 20 |
| 11774 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 |
| 12695 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 16166 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 19526 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 21741 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 24299 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 24339 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 17 of 20 |
| 24550 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 |
<!-- /table:evalFc_e0_pervoxel_1deg -->

<!-- table:evalFc_e0_pervoxel_2deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 121 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 |
| 1053 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 19 of 20 |
| 1078 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 14 of 20 |
| 2906 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 20 of 20 |
| 2910 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 6094 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 20 of 20 |
| 6674 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 20 of 20 |
| 7109 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 17 of 20 |
| 7577 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 16 of 20 |
| 8214 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 | 19 of 20 |
| 8814 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 11086 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 17 of 20 |
| 11774 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 18 of 20 |
| 12416 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 18 of 20 |
| 12695 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 18 of 20 |
| 16166 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 |
| 19526 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 21741 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 21831 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 |
| 22863 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 19 of 20 |
| 24299 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 |
| 24339 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 19 of 20 |
| 24550 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 15 of 20 |
<!-- /table:evalFc_e0_pervoxel_2deg -->

<!-- table:evalFc_e0_pervoxel_3deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 121 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 13 of 20 | 19 of 20 |
| 1053 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 13 of 20 |
| 1078 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 13 of 20 |
| 1260 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 16 of 20 |
| 2906 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 17 of 20 |
| 2910 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 14 of 20 |
| 6094 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 13 of 20 | 15 of 20 |
| 6674 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 18 of 20 |
| 7109 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 14 of 20 | 20 of 20 |
| 7577 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 17 of 20 |
| 8214 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 14 of 20 | 17 of 20 |
| 8814 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 14 of 20 | 13 of 20 |
| 11086 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 18 of 20 |
| 11774 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 19 of 20 |
| 12416 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 13 of 20 | 16 of 20 |
| 12695 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 15 of 20 |
| 16166 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 16 of 20 |
| 19526 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 13 of 20 | 19 of 20 |
| 21741 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 17 of 20 |
| 21831 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 12 of 20 | 13 of 20 |
| 22863 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 14 of 20 | 17 of 20 |
| 24299 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 15 of 20 |
| 24339 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 | 16 of 20 |
| 24550 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 17 of 20 |
<!-- /table:evalFc_e0_pervoxel_3deg -->

<!-- table:evalFc_e0_pervoxel_5deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 121 | 17 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 5 of 20 | 11 of 20 |
| 1053 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 6 of 20 | 6 of 20 |
| 1078 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 9 of 20 | 12 of 20 |
| 1260 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 10 of 20 | 11 of 20 |
| 2906 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 8 of 20 |
| 2910 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 9 of 20 | 12 of 20 |
| 6094 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 11 of 20 |
| 6674 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 12 of 20 | 13 of 20 |
| 7109 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 3 of 20 | 12 of 20 |
| 7577 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 3 of 20 | 14 of 20 |
| 8214 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 11 of 20 | 12 of 20 |
| 8814 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 9 of 20 | 11 of 20 |
| 11086 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 14 of 20 |
| 11774 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 11 of 20 |
| 12416 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 3 of 20 | 9 of 20 |
| 12695 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 10 of 20 | 14 of 20 |
| 16166 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 11 of 20 |
| 19526 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 11 of 20 | 8 of 20 |
| 21741 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 5 of 20 | 11 of 20 |
| 21831 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 6 of 20 | 10 of 20 |
| 22863 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 5 of 20 | 14 of 20 |
| 24299 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 10 of 20 |
| 24339 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 3 of 20 | 6 of 20 |
| 24550 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 10 of 20 |
<!-- /table:evalFc_e0_pervoxel_5deg -->

<!-- table:evalFc_e0_raw_1deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 13.0 | 9.00 | 268 | 342 |
| F4, 5-frame window | 2.00 | 0.00 | 13.0 | 9.00 | 274 | 391 |
| F5 expanding shells | 6.00 | 205 | 79.0 | 0.00 | 83.2 | 83.2 |
| F5 then F4 | 8.00 | 205 | 86.0 | 5.00 | 243 | 320 |
| Policy, 1-frame window | 2.00 | 0.00 | 13.0 | 9.00 | 268 | 342 |
| Policy, 5-frame window | 2.00 | 0.00 | 13.0 | 9.00 | 274 | 391 |
| CMA-ES, sigma0 1 deg | 281 | 0.00 | 0.00 | 0.00 | 281 | 358 |
| CMA-ES, sigma0 2 deg | 337 | 0.00 | 0.00 | 0.00 | 337 | 526 |
<!-- /table:evalFc_e0_raw_1deg -->

<!-- table:evalFc_e0_raw_2deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 15.0 | 11.0 | 304 | 385 |
| F4, 5-frame window | 2.00 | 0.00 | 18.5 | 11.0 | 304 | 429 |
| F5 expanding shells | 12.0 | 1324 | 175 | 0.00 | 219 | 279 |
| F5 then F4 | 14.0 | 1324 | 184 | 6.00 | 390 | 474 |
| Policy, 1-frame window | 2.00 | 0.00 | 15.0 | 11.0 | 304 | 385 |
| Policy, 5-frame window | 2.00 | 0.00 | 18.5 | 11.0 | 304 | 429 |
| CMA-ES, sigma0 1 deg | 323 | 0.00 | 0.00 | 0.00 | 323 | 435 |
| CMA-ES, sigma0 2 deg | 344 | 0.00 | 0.00 | 0.00 | 344 | 709 |
<!-- /table:evalFc_e0_raw_2deg -->

<!-- table:evalFc_e0_raw_3deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 18.0 | 12.0 | 331 | 404 |
| F4, 5-frame window | 2.00 | 0.00 | 22.0 | 12.0 | 331 | 447 |
| F5 expanding shells | 20.0 | 5593 | 303 | 0.00 | 516 | 516 |
| F5 then F4 | 22.0 | 5593 | 310 | 5.00 | 648 | 730 |
| Policy, 1-frame window | 2.00 | 0.00 | 18.0 | 12.0 | 331 | 412 |
| Policy, 5-frame window | 2.00 | 0.00 | 22.0 | 12.0 | 331 | 447 |
| CMA-ES, sigma0 1 deg | 376 | 0.00 | 0.00 | 0.00 | 376 | 911 |
| CMA-ES, sigma0 2 deg | 379 | 0.00 | 0.00 | 0.00 | 379 | 1109 |
<!-- /table:evalFc_e0_raw_3deg -->

<!-- table:evalFc_e0_raw_5deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 24.0 | 13.0 | 358 | 454 |
| F4, 5-frame window | 2.00 | 0.00 | 25.0 | 13.0 | 370 | 475 |
| F5 expanding shells | 32.0 | 21717 | 495 | 0.00 | 1364 | 1568 |
| F5 then F4 | 34.0 | 21717 | 504 | 6.00 | 1533 | 1709 |
| Policy, 1-frame window | 2.00 | 0.00 | 25.0 | 13.0 | 377 | 489 |
| Policy, 5-frame window | 2.00 | 0.00 | 25.0 | 13.0 | 370 | 475 |
| CMA-ES, sigma0 1 deg | 606 | 0.00 | 0.00 | 0.00 | 606 | 1591 |
| CMA-ES, sigma0 2 deg | 554 | 0.00 | 0.00 | 0.00 | 554 | 2220 |
<!-- /table:evalFc_e0_raw_5deg -->

<!-- table:evalFc_full_cheapest -->
| offset (deg) | share within 1 deg at least | fewest strict-equivalent evaluations (median) | fewest seconds per run (median) |
|---|---|---|---|
| 1 | 100% | F5 expanding shells (137) | F5 expanding shells (0.07) |
| 1 | 99% | F5 expanding shells (137) | F5 expanding shells (0.07) |
| 1 | 95% | F5 expanding shells (137) | F5 expanding shells (0.07) |
| 1 | 90% | F5 expanding shells (137) | F5 expanding shells (0.07) |
| 2 | 100% | Policy, 1-frame window (379) | Policy, 1-frame window (0.14) |
| 2 | 99% | Policy, 1-frame window (379) | Policy, 1-frame window (0.14) |
| 2 | 95% | CMA-ES, sigma0 1 deg (288) | CMA-ES, sigma0 1 deg (0.11) |
| 2 | 90% | CMA-ES, sigma0 1 deg (288) | CMA-ES, sigma0 1 deg (0.11) |
| 3 | 100% | Policy, 5-frame window (410) | Policy, 1-frame window (0.16) |
| 3 | 99% | Policy, 5-frame window (410) | Policy, 1-frame window (0.16) |
| 3 | 95% | Policy, 5-frame window (410) | Policy, 1-frame window (0.16) |
| 3 | 90% | Policy, 5-frame window (410) | Policy, 1-frame window (0.16) |
| 5 | 100% | Policy, 5-frame window (530) | Policy, 5-frame window (0.18) |
| 5 | 99% | Policy, 5-frame window (530) | Policy, 5-frame window (0.18) |
| 5 | 95% | Policy, 5-frame window (530) | Policy, 5-frame window (0.18) |
| 5 | 90% | Policy, 5-frame window (530) | Policy, 5-frame window (0.18) |
<!-- /table:evalFc_full_cheapest -->

<!-- table:evalFc_full_contention -->
| offset (deg) | option | single-worker, median (s) | contended, same runs, median (s) | ratio of medians |
|---|---|---|---|---|
| 1 | F4, 1-frame window | 0.12 | 0.14 | 1.17 |
| 1 | F4, 5-frame window | 0.12 | 0.14 | 1.13 |
| 1 | F5 expanding shells | 0.07 | 0.08 | 1.11 |
| 1 | F5 then F4 | 0.17 | 0.19 | 1.10 |
| 1 | Policy, 1-frame window | 0.12 | 0.14 | 1.17 |
| 1 | Policy, 5-frame window | 0.12 | 0.14 | 1.13 |
| 1 | CMA-ES, sigma0 1 deg | 0.11 | 0.13 | 1.20 |
| 1 | CMA-ES, sigma0 2 deg | 0.13 | 0.15 | 1.11 |
| 2 | F4, 1-frame window | 0.14 | 0.16 | 1.12 |
| 2 | F4, 5-frame window | 0.15 | 0.16 | 1.13 |
| 2 | F5 expanding shells | 0.27 | 0.30 | 1.10 |
| 2 | F5 then F4 | 0.35 | 0.38 | 1.09 |
| 2 | Policy, 1-frame window | 0.14 | 0.16 | 1.12 |
| 2 | Policy, 5-frame window | 0.15 | 0.16 | 1.13 |
| 2 | CMA-ES, sigma0 1 deg | 0.11 | 0.13 | 1.13 |
| 2 | CMA-ES, sigma0 2 deg | 0.13 | 0.15 | 1.17 |
| 3 | F4, 1-frame window | 0.16 | 0.18 | 1.12 |
| 3 | F4, 5-frame window | 0.16 | 0.18 | 1.12 |
| 3 | F5 expanding shells | 0.46 | 0.50 | 1.08 |
| 3 | F5 then F4 | 0.54 | 0.59 | 1.08 |
| 3 | Policy, 1-frame window | 0.16 | 0.18 | 1.12 |
| 3 | Policy, 5-frame window | 0.16 | 0.18 | 1.12 |
| 3 | CMA-ES, sigma0 1 deg | 0.15 | 0.17 | 1.10 |
| 3 | CMA-ES, sigma0 2 deg | 0.18 | 0.20 | 1.15 |
| 5 | F4, 1-frame window | 0.17 | 0.19 | 1.11 |
| 5 | F4, 5-frame window | 0.18 | 0.19 | 1.07 |
| 5 | F5 expanding shells | 1.05 | 1.18 | 1.12 |
| 5 | F5 then F4 | 1.14 | 1.28 | 1.12 |
| 5 | Policy, 1-frame window | 0.19 | 0.20 | 1.08 |
| 5 | Policy, 5-frame window | 0.18 | 0.20 | 1.11 |
| 5 | CMA-ES, sigma0 1 deg | 0.20 | 0.22 | 1.13 |
| 5 | CMA-ES, sigma0 2 deg | 0.16 | 0.18 | 1.13 |
<!-- /table:evalFc_full_contention -->

<!-- table:evalFc_full_main_1deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 414 of 480, 86.2% (82.9-89.0) | n/a | 317 | 427 | 0.12 | 0.17 |
| F4, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 372 of 480, 77.5% (73.6-81.0) | n/a | 317 | 413 | 0.12 | 0.15 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 137 | 138 | 0.07 | 0.08 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 411 of 480, 85.6% (82.2-88.5) | n/a | 392 | 483 | 0.17 | 0.19 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 414 of 480, 86.2% (82.9-89.0) | 0 of 480 | 317 | 427 | 0.12 | 0.17 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 372 of 480, 77.5% (73.6-81.0) | 0 of 480 | 317 | 413 | 0.12 | 0.15 |
| CMA-ES, sigma0 1 deg | 480 of 480, 100.0% (99.2-100.0) | 476 of 480, 99.2% (97.9-99.7) | n/a | 267 | 337 | 0.11 | 0.18 |
| CMA-ES, sigma0 2 deg | 466 of 480, 97.1% (95.2-98.3) | 463 of 480, 96.5% (94.4-97.8) | n/a | 309 | 414 | 0.13 | 0.15 |
<!-- /table:evalFc_full_main_1deg -->

<!-- table:evalFc_full_main_2deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 479 of 480, 99.8% (98.8-100.0) | 414 of 480, 86.2% (82.9-89.0) | n/a | 379 | 482 | 0.14 | 0.17 |
| F4, 5-frame window | 479 of 480, 99.8% (98.8-100.0) | 375 of 480, 78.1% (74.2-81.6) | n/a | 379 | 471 | 0.15 | 0.17 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 437 | 548 | 0.27 | 0.28 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 410 of 480, 85.4% (82.0-88.3) | n/a | 738 | 834 | 0.35 | 0.37 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 415 of 480, 86.5% (83.1-89.2) | 1 of 480 | 379 | 483 | 0.14 | 0.17 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 376 of 480, 78.3% (74.4-81.8) | 1 of 480 | 379 | 471 | 0.15 | 0.17 |
| CMA-ES, sigma0 1 deg | 468 of 480, 97.5% (95.7-98.6) | 462 of 480, 96.2% (94.2-97.6) | n/a | 288 | 372 | 0.11 | 0.12 |
| CMA-ES, sigma0 2 deg | 446 of 480, 92.9% (90.3-94.9) | 440 of 480, 91.7% (88.9-93.8) | n/a | 323 | 450 | 0.13 | 0.23 |
<!-- /table:evalFc_full_main_2deg -->

<!-- table:evalFc_full_main_3deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 449 of 480, 93.5% (91.0-95.4) | 387 of 480, 80.6% (76.9-83.9) | n/a | 410 | 536 | 0.16 | 0.18 |
| F4, 5-frame window | 466 of 480, 97.1% (95.2-98.3) | 363 of 480, 75.6% (71.6-79.3) | n/a | 410 | 530 | 0.16 | 0.18 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 931 | 931 | 0.46 | 0.47 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 411 of 480, 85.6% (82.2-88.5) | n/a | 1146 | 1249 | 0.54 | 0.57 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 415 of 480, 86.5% (83.1-89.2) | 31 of 480 | 417 | 595 | 0.16 | 0.18 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 373 of 480, 77.7% (73.8-81.2) | 14 of 480 | 410 | 537 | 0.16 | 0.18 |
| CMA-ES, sigma0 1 deg | 395 of 480, 82.3% (78.6-85.4) | 388 of 480, 80.8% (77.1-84.1) | n/a | 337 | 457 | 0.15 | 0.21 |
| CMA-ES, sigma0 2 deg | 415 of 480, 86.5% (83.1-89.2) | 407 of 480, 84.8% (81.3-87.7) | n/a | 351 | 540 | 0.18 | 0.23 |
<!-- /table:evalFc_full_main_3deg -->

<!-- table:evalFc_full_main_5deg -->
| option | within 1 deg | within 0.05 deg | fell back | strict-equivalent, median | strict-equivalent, 90th pct | seconds per run, median | seconds per run, 90th pct |
|---|---|---|---|---|---|---|---|
| F4, 1-frame window | 315 of 480, 65.6% (61.3-69.7) | 279 of 480, 58.1% (53.7-62.5) | n/a | 461 | 657 | 0.17 | 0.20 |
| F4, 5-frame window | 344 of 480, 71.7% (67.5-75.5) | 274 of 480, 57.1% (52.6-61.4) | n/a | 472 | 668 | 0.18 | 0.26 |
| F5 expanding shells | 480 of 480, 100.0% (99.2-100.0) | 0 of 480, 0.0% (0.0-0.8) | n/a | 2026 | 2265 | 1.05 | 1.09 |
| F5 then F4 | 480 of 480, 100.0% (99.2-100.0) | 409 of 480, 85.2% (81.8-88.1) | n/a | 2353 | 2547 | 1.14 | 1.17 |
| Policy, 1-frame window | 480 of 480, 100.0% (99.2-100.0) | 410 of 480, 85.4% (82.0-88.3) | 165 of 480 | 562 | 2841 | 0.19 | 1.28 |
| Policy, 5-frame window | 480 of 480, 100.0% (99.2-100.0) | 368 of 480, 76.7% (72.7-80.2) | 136 of 480 | 530 | 2891 | 0.18 | 0.36 |
| CMA-ES, sigma0 1 deg | 165 of 480, 34.4% (30.3-38.7) | 162 of 480, 33.8% (29.7-38.1) | n/a | 421 | 659 | 0.20 | 0.40 |
| CMA-ES, sigma0 2 deg | 272 of 480, 56.7% (52.2-61.0) | 261 of 480, 54.4% (49.9-58.8) | n/a | 428 | 785 | 0.16 | 0.29 |
<!-- /table:evalFc_full_main_5deg -->

<!-- table:evalFc_full_paired -->
| offset (deg) | comparison | fixed | broken | both wrong | both right | voxels with a fixed run | voxels with a broken run |
|---|---|---|---|---|---|---|---|
| 1 | Policy 1 against F4 (1-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 1 | Policy 5 against F4 (5-frame) | 0 | 0 | 0 | 480 | 0 | 0 |
| 1 | F4 5-frame against F4 1-frame | 0 | 0 | 0 | 480 | 0 | 0 |
| 1 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
| 2 | Policy 1 against F4 (1-frame) | 1 | 0 | 0 | 479 | 1 | 0 |
| 2 | Policy 5 against F4 (5-frame) | 1 | 0 | 0 | 479 | 1 | 0 |
| 2 | F4 5-frame against F4 1-frame | 0 | 0 | 1 | 479 | 0 | 0 |
| 2 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
| 3 | Policy 1 against F4 (1-frame) | 31 | 0 | 0 | 449 | 17 | 0 |
| 3 | Policy 5 against F4 (5-frame) | 14 | 0 | 0 | 466 | 10 | 0 |
| 3 | F4 5-frame against F4 1-frame | 17 | 0 | 14 | 449 | 14 | 0 |
| 3 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
| 5 | Policy 1 against F4 (1-frame) | 165 | 0 | 0 | 315 | 24 | 0 |
| 5 | Policy 5 against F4 (5-frame) | 136 | 0 | 0 | 344 | 24 | 0 |
| 5 | F4 5-frame against F4 1-frame | 33 | 4 | 132 | 311 | 15 | 3 |
| 5 | Policy 5 against policy 1 | 0 | 0 | 0 | 480 | 0 | 0 |
<!-- /table:evalFc_full_paired -->

<!-- table:evalFc_full_pervoxel_1deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 704 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 3156 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 9861 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 9864 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 11924 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 11927 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 12263 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 13333 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 13455 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 14489 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 14771 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 18535 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 19577 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 22800 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
<!-- /table:evalFc_full_pervoxel_1deg -->

<!-- table:evalFc_full_pervoxel_2deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 704 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 1729 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 18 of 20 |
| 3156 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 18 of 20 |
| 3188 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 3633 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 3634 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 20 of 20 |
| 9861 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 |
| 9864 | 19 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
| 11924 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 18 of 20 |
| 11927 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 12263 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 13333 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 20 of 20 |
| 13455 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 |
| 14489 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 |
| 14492 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 17 of 20 |
| 14771 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 19 of 20 |
| 17490 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 18 of 20 |
| 18535 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 20 of 20 |
| 19572 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 20 of 20 |
| 23300 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 |
<!-- /table:evalFc_full_pervoxel_2deg -->

<!-- table:evalFc_full_pervoxel_3deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 704 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 17 of 20 |
| 1729 | 19 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 18 of 20 |
| 3156 | 17 of 20 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 | 15 of 20 |
| 3188 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 16 of 20 |
| 3285 | 17 of 20 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 19 of 20 |
| 3633 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 20 of 20 |
| 3634 | 17 of 20 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 13 of 20 | 16 of 20 |
| 9861 | 19 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 19 of 20 |
| 9864 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 19 of 20 |
| 11924 | 16 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 12 of 20 | 15 of 20 |
| 11927 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 | 16 of 20 |
| 12260 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 18 of 20 | 20 of 20 |
| 12263 | 16 of 20 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 17 of 20 |
| 13333 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 13 of 20 | 15 of 20 |
| 13455 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 14 of 20 | 18 of 20 |
| 14489 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 15 of 20 | 13 of 20 |
| 14492 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 14 of 20 |
| 14771 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 18 of 20 |
| 17490 | 18 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 | 19 of 20 |
| 18535 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 19 of 20 |
| 19572 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 17 of 20 | 19 of 20 |
| 19577 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 18 of 20 |
| 22800 | 19 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 19 of 20 | 19 of 20 |
| 23300 | 18 of 20 | 19 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 16 of 20 | 16 of 20 |
<!-- /table:evalFc_full_pervoxel_3deg -->

<!-- table:evalFc_full_pervoxel_5deg -->
| voxel | F4 | F4w | F5 | F54 | P1 | P5 | C1 | C2 |
|---|---|---|---|---|---|---|---|---|
| 704 | 15 of 20 | 17 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 3 of 20 | 13 of 20 |
| 1729 | 16 of 20 | 16 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 4 of 20 | 10 of 20 |
| 3156 | 11 of 20 | 12 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 4 of 20 | 14 of 20 |
| 3188 | 15 of 20 | 15 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 4 of 20 | 12 of 20 |
| 3285 | 9 of 20 | 9 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 9 of 20 | 15 of 20 |
| 3633 | 12 of 20 | 12 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 5 of 20 |
| 3634 | 11 of 20 | 15 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 10 of 20 | 12 of 20 |
| 9861 | 10 of 20 | 11 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 5 of 20 | 13 of 20 |
| 9864 | 15 of 20 | 15 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 11 of 20 |
| 11924 | 10 of 20 | 13 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 6 of 20 | 12 of 20 |
| 11927 | 15 of 20 | 14 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 7 of 20 |
| 12260 | 16 of 20 | 16 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 13 of 20 |
| 12263 | 5 of 20 | 7 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 6 of 20 | 11 of 20 |
| 13333 | 14 of 20 | 16 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 6 of 20 | 15 of 20 |
| 13455 | 16 of 20 | 16 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 13 of 20 |
| 14489 | 14 of 20 | 17 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 9 of 20 | 10 of 20 |
| 14492 | 16 of 20 | 17 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 11 of 20 |
| 14771 | 15 of 20 | 18 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 9 of 20 |
| 17490 | 12 of 20 | 12 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 11 of 20 |
| 18535 | 17 of 20 | 17 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 9 of 20 | 13 of 20 |
| 19572 | 8 of 20 | 13 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 14 of 20 |
| 19577 | 13 of 20 | 13 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 7 of 20 | 14 of 20 |
| 22800 | 14 of 20 | 16 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 7 of 20 |
| 23300 | 16 of 20 | 17 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 8 of 20 | 7 of 20 |
<!-- /table:evalFc_full_pervoxel_5deg -->

<!-- table:evalFc_full_raw_1deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 17.0 | 10.0 | 317 | 427 |
| F4, 5-frame window | 2.00 | 0.00 | 16.0 | 10.0 | 317 | 413 |
| F5 expanding shells | 4.00 | 74.0 | 47.0 | 0.00 | 137 | 138 |
| F5 then F4 | 6.00 | 74.0 | 59.5 | 7.00 | 392 | 483 |
| Policy, 1-frame window | 2.00 | 0.00 | 17.0 | 10.0 | 317 | 427 |
| Policy, 5-frame window | 2.00 | 0.00 | 16.0 | 10.0 | 317 | 413 |
| CMA-ES, sigma0 1 deg | 267 | 0.00 | 0.00 | 0.00 | 267 | 337 |
| CMA-ES, sigma0 2 deg | 309 | 0.00 | 0.00 | 0.00 | 309 | 414 |
<!-- /table:evalFc_full_raw_1deg -->

<!-- table:evalFc_full_raw_2deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 20.0 | 11.0 | 379 | 482 |
| F4, 5-frame window | 2.00 | 0.00 | 17.0 | 11.0 | 379 | 471 |
| F5 expanding shells | 10.0 | 801 | 143 | 0.00 | 437 | 548 |
| F5 then F4 | 12.0 | 801 | 162 | 8.00 | 738 | 834 |
| Policy, 1-frame window | 2.00 | 0.00 | 20.0 | 11.0 | 379 | 483 |
| Policy, 5-frame window | 2.00 | 0.00 | 17.0 | 11.0 | 379 | 471 |
| CMA-ES, sigma0 1 deg | 288 | 0.00 | 0.00 | 0.00 | 288 | 372 |
| CMA-ES, sigma0 2 deg | 323 | 0.00 | 0.00 | 0.00 | 323 | 450 |
<!-- /table:evalFc_full_raw_2deg -->

<!-- table:evalFc_full_raw_3deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 20.0 | 12.0 | 410 | 536 |
| F4, 5-frame window | 2.00 | 0.00 | 19.0 | 13.0 | 410 | 530 |
| F5 expanding shells | 18.0 | 4141 | 271 | 0.00 | 931 | 931 |
| F5 then F4 | 20.0 | 4141 | 279 | 8.00 | 1146 | 1249 |
| Policy, 1-frame window | 2.00 | 0.00 | 22.0 | 13.0 | 417 | 595 |
| Policy, 5-frame window | 2.00 | 0.00 | 19.0 | 13.0 | 410 | 537 |
| CMA-ES, sigma0 1 deg | 337 | 0.00 | 0.00 | 0.00 | 337 | 457 |
| CMA-ES, sigma0 2 deg | 351 | 0.00 | 0.00 | 0.00 | 351 | 540 |
<!-- /table:evalFc_full_raw_3deg -->

<!-- table:evalFc_full_raw_5deg -->
| option | strict | screen | smooth value | smooth Jacobian | strict-equivalent, median | strict-equivalent, 90th pct |
|---|---|---|---|---|---|---|
| F4, 1-frame window | 2.00 | 0.00 | 23.0 | 14.0 | 461 | 657 |
| F4, 5-frame window | 2.00 | 0.00 | 22.0 | 14.0 | 472 | 668 |
| F5 expanding shells | 30.0 | 18002 | 463 | 0.00 | 2026 | 2265 |
| F5 then F4 | 32.0 | 18002 | 483 | 8.00 | 2353 | 2547 |
| Policy, 1-frame window | 2.00 | 0.00 | 28.0 | 16.0 | 562 | 2841 |
| Policy, 5-frame window | 2.00 | 0.00 | 26.0 | 16.0 | 530 | 2891 |
| CMA-ES, sigma0 1 deg | 421 | 0.00 | 0.00 | 0.00 | 421 | 659 |
| CMA-ES, sigma0 2 deg | 428 | 0.00 | 0.00 | 0.00 | 428 | 785 |
<!-- /table:evalFc_full_raw_5deg -->
