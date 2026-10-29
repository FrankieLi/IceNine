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
