<!-- table:evalE_e2e -->
| voxel set | rule | arm (SmoothCostRank) | n | wrong, baseline | wrong, arm | fixed | broken | both wrong | broken: traced to arm / draw / untraced | fixed: traced to arm / draw | median error of voxels right in both (deg), baseline / arm | strict evaluations per voxel, baseline | strict evaluations per voxel, arm | smooth evaluations per voxel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E0 200 | lenient | E3+E4,tau=0.25 | 200 | 41 (15.5-26.6%) | 7 (1.7-7.0%) | 34 | 0 | 7 | 0 / 0 / 0 | 17 / 17 | 0.0130 / 0.0054 | 48700 | 49215 | 252 |
| E0 200 | C++ | E3+E4,tau=0.25 | 200 | 56 (22.2-34.6%) | 9 (2.4-8.3%) | 47 | 0 | 9 | 0 / 0 / 0 | 21 / 26 | 0.0223 / 0.0082 | 47230 | 47551 | 172 |
| six full-sample | lenient | E3+E4,tau=0.25 | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0181 / 0.0163 | 332046 | 334486 | 15030 |
| six full-sample | C++ | E3+E4,tau=0.25 | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0286 / 0.0291 | 228146 | 235524 | 11084 |
| E0 200 | lenient | E1+E2+E3+E4 | 200 | 41 (15.5-26.6%) | 5 (1.1-5.7%) | 36 | 0 | 5 | 0 / 0 / 0 | 15 / 21 | 0.0130 / 0.0050 | 48700 | 49653 | 263 |
| E0 200 | C++ | E1+E2+E3+E4 | 200 | 56 (22.2-34.6%) | 9 (2.4-8.3%) | 47 | 0 | 9 | 0 / 0 / 0 | 21 / 26 | 0.0223 / 0.0112 | 47230 | 47767 | 177 |
| six full-sample | lenient | E1+E2+E3+E4 | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0181 / 0.0284 | 332046 | 348018 | 15755 |
<!-- /table:evalE_e2e -->

<!-- table:evalE_landL0_basin_e0 -->
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays monotone to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays monotone to 5 deg, eps 0.05 |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.113 | 5.000 | 5.000 | 598 of 600 | 5.000 | 600 of 600 |
| strict cost, lenient rule | 0.113 | 5.000 | 5.000 | 599 of 600 | 5.000 | 600 of 600 |
| E1, level-0 width (16.5 px) | 0.843 | 5.000 | 5.000 | 566 of 600 | 5.000 | 600 of 600 |
| E2, level-0 reach (1.67 frames) | 0.142 | 5.000 | 5.000 | 599 of 600 | 5.000 | 600 of 600 |
| E4, level-0 tau (0.42 deg) | 1.449 | 5.000 | 5.000 | 548 of 600 | 5.000 | 600 of 600 |
| E1 + E2, level-0 widths | 1.975 | 5.000 | 5.000 | 528 of 600 | 5.000 | 600 of 600 |
| E1 + E2 + E3, level-0 widths | 1.975 | 5.000 | 5.000 | 592 of 600 | 5.000 | 600 of 600 |
| E2 + E3, level-0 reach | 0.142 | 5.000 | 5.000 | 600 of 600 | 5.000 | 600 of 600 |
| E3 + E4, level-0 tau | 0.454 | 5.000 | 5.000 | 598 of 600 | 5.000 | 600 of 600 |
| E1 + E2 + E3 + E4, level-0 widths | 1.828 | 5.000 | 5.000 | 592 of 600 | 5.000 | 600 of 600 |
<!-- /table:evalE_landL0_basin_e0 -->

<!-- table:evalE_landL0_basin_full -->
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays monotone to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays monotone to 5 deg, eps 0.05 |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.312 | 3.693 | 1.729 | 45 of 120 | 5.000 | 120 of 120 |
| strict cost, lenient rule | 0.333 | 2.783 | 1.971 | 15 of 120 | 5.000 | 120 of 120 |
| E1, level-0 width (16.5 px) | 1.795 | 1.006 | 0.766 | 0 of 120 | 4.412 | 55 of 120 |
| E2, level-0 reach (1.67 frames) | 0.517 | 2.002 | 1.375 | 4 of 120 | 5.000 | 111 of 120 |
| E4, level-0 tau (0.42 deg) | 1.434 | 5.000 | 1.519 | 77 of 120 | 5.000 | 120 of 120 |
| E1 + E2, level-0 widths | beyond 5 | 3.484 | 0.971 | 25 of 120 | 5.000 | 115 of 120 |
| E1 + E2 + E3, level-0 widths | beyond 5 | 4.655 | 3.690 | 49 of 120 | 5.000 | 112 of 120 |
| E2 + E3, level-0 reach | 0.520 | 2.121 | 1.718 | 3 of 120 | 5.000 | 103 of 120 |
| E3 + E4, level-0 tau | 0.610 | 5.000 | 4.011 | 67 of 120 | 5.000 | 120 of 120 |
| E1 + E2 + E3 + E4, level-0 widths | 2.602 | 5.000 | 4.820 | 87 of 120 | 5.000 | 120 of 120 |
<!-- /table:evalE_landL0_basin_full -->

<!-- table:evalE_landL0_minima_e0 -->
| cost | rays with a minimum of depth 0.01, dense, 0 to 0.5 deg | rays with a minimum of depth 0.05, dense, 0 to 0.5 deg | rays where the truth is not the lowest point | lowest point when not the truth: median deg (max) | ... and its depth below the truth (median) |
|---|---|---|---|---|---|
| strict cost, C++ rule | 12 of 600 (1.1-3.5%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0005 (0.0010) | 0.0021 |
| strict cost, lenient rule | 13 of 600 (1.3-3.7%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0005 (0.0010) | 0.0021 |
| E1, level-0 width (16.5 px) | 61 of 600 (8.0-12.8%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E2, level-0 reach (1.67 frames) | 17 of 600 (1.8-4.5%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E4, level-0 tau (0.42 deg) | 70 of 600 (9.3-14.5%) | 0 of 600 (0.0-0.6%) | 292 of 600 | 0.0085 (0.1730) | 0.0000 |
| E1 + E2, level-0 widths | 83 of 600 (11.3-16.8%) | 0 of 600 (0.0-0.6%) | 49 of 600 | 0.1370 (0.4475) | 0.0119 |
| E1 + E2 + E3, level-0 widths | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 30 of 600 | 0.0235 (0.0450) | 0.0000 |
| E2 + E3, level-0 reach | 0 of 600 (0.0-0.6%) | 0 of 600 (0.0-0.6%) | 6 of 600 | 0.0008 (0.0030) | 0.0000 |
| E3 + E4, level-0 tau | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 124 of 600 | 0.0010 (0.0060) | 0.0000 |
| E1 + E2 + E3 + E4, level-0 widths | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 281 of 600 | 0.0080 (0.0330) | 0.0000 |
<!-- /table:evalE_landL0_minima_e0 -->

<!-- table:evalE_landL0_minima_full -->
| cost | rays with a minimum of depth 0.01, dense, 0 to 5 deg | rays with a minimum of depth 0.05, dense, 0 to 5 deg | rays with a minimum of depth 0.01, 0 to 0.5 deg | rays with a minimum of depth 0.05, 0 to 0.5 deg | rays where the truth is not the lowest point | lowest point when not the truth: median deg (max) | ... and its depth below the truth (median) |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 72 of 120 (51.1-68.3%) | 0 of 120 (0.0-3.1%) | 6 of 120 (2.3-10.5%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| strict cost, lenient rule | 97 of 120 (72.9-86.9%) | 0 of 120 (0.0-3.1%) | 5 of 120 (1.8-9.4%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E1, level-0 width (16.5 px) | 120 of 120 (96.9-100.0%) | 46 of 120 (30.1-47.3%) | 5 of 120 (1.8-9.4%) | 0 of 120 (0.0-3.1%) | 36 of 120 | 0.0203 (0.0460) | 0.0000 |
| E2, level-0 reach (1.67 frames) | 113 of 120 (88.4-97.1%) | 6 of 120 (2.3-10.5%) | 18 of 120 (9.7-22.5%) | 0 of 120 (0.0-3.1%) | 2 of 120 | 0.2083 (0.2230) | 0.0057 |
| E4, level-0 tau (0.42 deg) | 39 of 120 (24.8-41.3%) | 0 of 120 (0.0-3.1%) | 24 of 120 (13.8-28.0%) | 0 of 120 (0.0-3.1%) | 55 of 120 | 0.0235 (0.2230) | 0.0001 |
| E1 + E2, level-0 widths | 74 of 120 (52.7-69.9%) | 0 of 120 (0.0-3.1%) | 16 of 120 (8.4-20.6%) | 0 of 120 (0.0-3.1%) | 35 of 120 | 0.0960 (0.2950) | 0.0019 |
| E1 + E2 + E3, level-0 widths | 34 of 120 (21.0-37.0%) | 1 of 120 (0.1-4.6%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E2 + E3, level-0 reach | 115 of 120 (90.6-98.2%) | 10 of 120 (4.6-14.7%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E3 + E4, level-0 tau | 46 of 120 (30.1-47.3%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 65 of 120 | 0.0145 (0.0575) | 0.0000 |
| E1 + E2 + E3 + E4, level-0 widths | 16 of 120 (8.4-20.6%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 72 of 120 | 0.0255 (0.0765) | 0.0000 |
<!-- /table:evalE_landL0_minima_full -->

<!-- table:evalE_landL0_steps_e0 -->
| cost | cost at the truth (median) | steps per deg, under 0.05 deg | steps per deg, 0.05 to 0.5 deg | steps per deg, 0.5 to 5 deg (4 windows of 0.05 deg) | variation in steps, under 0.05 deg | variation in steps, 0.05 to 0.5 deg | variation in steps, 0.5 to 5 deg (4 windows of 0.05 deg) |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.0606 | 330.8 | 284.7 | 25.8 | 40.7% | 61.9% | 86.2% |
| strict cost, lenient rule | 0.0606 | 318.5 | 349.5 | 50.2 | 38.0% | 55.5% | 84.2% |
| E1, level-0 width (16.5 px) | 0.0606 | 42.7 | 42.3 | 14.2 | 96.1% | 93.1% | 87.2% |
| E2, level-0 reach (1.67 frames) | 0.0616 | 217.8 | 235.2 | 39.0 | 21.0% | 38.0% | 73.4% |
| E4, level-0 tau (0.42 deg) | 0.0018 | 1.8 | 4.3 | 11.9 | 39.1% | 9.0% | 13.2% |
| E1 + E2, level-0 widths | 0.0606 | 2.2 | 1.8 | 1.4 | 72.1% | 14.5% | 4.6% |
| E1 + E2 + E3, level-0 widths | 0.0000 | 0.0 | 0.0 | 0.0 | 0.0% | 0.6% | 0.1% |
| E2 + E3, level-0 reach | 0.0000 | 238.0 | 241.2 | 38.1 | 23.5% | 40.0% | 73.9% |
| E3 + E4, level-0 tau | 0.0008 | 250.7 | 230.3 | 35.2 | 31.3% | 43.2% | 32.9% |
| E1 + E2 + E3 + E4, level-0 widths | 0.0008 | 0.1 | 2.8 | 9.1 | 0.9% | 2.4% | 5.5% |
<!-- /table:evalE_landL0_steps_e0 -->

<!-- table:evalE_landL0_steps_full -->
| cost | cost at the truth (median) | steps per deg, under 0.05 deg | steps per deg, 0.05 to 0.5 deg | steps per deg, 0.5 to 5 deg | variation in steps, under 0.05 deg | variation in steps, 0.05 to 0.5 deg | variation in steps, 0.5 to 5 deg |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.0703 | 113.8 | 261.5 | 25.0 | 88.2% | 84.2% | 94.0% |
| strict cost, lenient rule | 0.0703 | 116.3 | 323.6 | 108.2 | 87.5% | 78.4% | 89.6% |
| E1, level-0 width (16.5 px) | 0.0671 | 49.3 | 49.7 | 49.0 | 97.3% | 94.6% | 91.4% |
| E2, level-0 reach (1.67 frames) | 0.0717 | 55.2 | 239.3 | 113.8 | 56.1% | 55.6% | 69.0% |
| E4, level-0 tau (0.42 deg) | 0.2820 | 6.8 | 5.0 | 15.2 | 1.9% | 11.8% | 18.1% |
| E1 + E2, level-0 widths | 0.0682 | 0.0 | 2.3 | 2.6 | 0.0% | 18.1% | 12.7% |
| E1 + E2 + E3, level-0 widths | 0.0000 | 0.0 | 0.0 | 0.0 | 0.0% | 0.0% | 0.3% |
| E2 + E3, level-0 reach | 0.0000 | 56.0 | 241.2 | 112.2 | 57.9% | 55.9% | 69.1% |
| E3 + E4, level-0 tau | 0.1407 | 96.3 | 240.4 | 56.4 | 85.7% | 66.6% | 57.1% |
| E1 + E2 + E3 + E4, level-0 widths | 0.1407 | 6.3 | 3.4 | 11.7 | 1.4% | 2.5% | 8.7% |
<!-- /table:evalE_landL0_steps_full -->

<!-- table:evalE_landL0_types_e0 -->
| cost | isolated steps | at a frame event | at a pixel event | at a presence event | at a detector event | at no strict event |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 89902 | 10.4% | 89.1% | 0.4% | 0.1% | 0.0% |
| strict cost, lenient rule | 109934 | 10.0% | 89.6% | 0.3% | 0.1% | 0.0% |
| E1, level-0 width (16.5 px) | 14414 | 93.6% | 1.5% | 3.5% | 1.2% | 0.2% |
| E2, level-0 reach (1.67 frames) | 74716 | 4.1% | 90.5% | 0.5% | 0.2% | 4.7% |
| E4, level-0 tau (0.42 deg) | 2646 | 74.9% | 6.9% | 17.4% | 0.1% | 0.6% |
| E1 + E2, level-0 widths | 709 | 2.4% | 3.0% | 67.6% | 27.1% | 0.0% |
| E1 + E2 + E3, level-0 widths | 7 | 0.0% | 0.0% | 100.0% | 0.0% | 0.0% |
| E2 + E3, level-0 reach | 76840 | 4.2% | 91.0% | 0.1% | 0.1% | 4.6% |
| E3 + E4, level-0 tau | 73931 | 15.2% | 80.7% | 0.1% | 0.1% | 3.9% |
| E1 + E2 + E3 + E4, level-0 widths | 1863 | 97.4% | 1.7% | 0.6% | 0.0% | 0.3% |
<!-- /table:evalE_landL0_types_e0 -->

<!-- table:evalE_landL0_types_full -->
| cost | isolated steps | at a frame event | at a pixel event | at a presence event | at a detector event | at no strict event |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 28290 | 15.2% | 83.1% | 1.6% | 0.1% | 0.0% |
| strict cost, lenient rule | 76615 | 10.6% | 88.5% | 0.8% | 0.1% | 0.0% |
| E1, level-0 width (16.5 px) | 29433 | 92.5% | 2.4% | 2.5% | 2.4% | 0.2% |
| E2, level-0 reach (1.67 frames) | 74698 | 4.0% | 91.4% | 0.8% | 0.4% | 3.3% |
| E4, level-0 tau (0.42 deg) | 8497 | 87.9% | 4.2% | 7.3% | 0.2% | 0.4% |
| E1 + E2, level-0 widths | 1531 | 2.4% | 9.0% | 42.6% | 45.2% | 0.8% |
| E1 + E2 + E3, level-0 widths | 14 | 0.0% | 7.1% | 85.7% | 7.1% | 0.0% |
| E2 + E3, level-0 reach | 73943 | 4.1% | 92.4% | 0.1% | 0.1% | 3.3% |
| E3 + E4, level-0 tau | 44000 | 30.1% | 67.3% | 0.2% | 0.1% | 2.2% |
| E1 + E2 + E3 + E4, level-0 widths | 6552 | 98.0% | 1.2% | 0.5% | 0.0% | 0.3% |
<!-- /table:evalE_landL0_types_full -->

<!-- table:evalE_land_basin_e0 -->
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays monotone to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays monotone to 5 deg, eps 0.05 |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.113 | 5.000 | 5.000 | 598 of 600 | 5.000 | 600 of 600 |
| strict cost, lenient rule | 0.113 | 5.000 | 5.000 | 599 of 600 | 5.000 | 600 of 600 |
| S: sampled binary (control) | 0.122 | 5.000 | 5.000 | 600 of 600 | 5.000 | 600 of 600 |
| E1: distance kernel 3 px | 0.530 | 5.000 | 5.000 | 577 of 600 | 5.000 | 600 of 600 |
| E1: distance kernel 12 px | 0.812 | 5.000 | 5.000 | 567 of 600 | 5.000 | 600 of 600 |
| E2: omega blending, 1 frame | 0.132 | 5.000 | 5.000 | 599 of 600 | 5.000 | 600 of 600 |
| E3: soft peak weights | 0.122 | 5.000 | 5.000 | 600 of 600 | 5.000 | 600 of 600 |
| E4: ray residual, tau 0.05 deg | 0.454 | 5.000 | 5.000 | 588 of 600 | 5.000 | 600 of 600 |
| E4: ray residual, tau 0.25 deg | 1.149 | 5.000 | 5.000 | 559 of 600 | 5.000 | 600 of 600 |
| E1 + E2 | 0.843 | 5.000 | 5.000 | 558 of 600 | 5.000 | 600 of 600 |
| E1 + E2 + E3 | 0.843 | 5.000 | 5.000 | 597 of 600 | 5.000 | 600 of 600 |
| E2 + E3 | 0.142 | 5.000 | 5.000 | 600 of 600 | 5.000 | 600 of 600 |
| E3 + E4 (0.25 deg) | 0.389 | 5.000 | 5.000 | 600 of 600 | 5.000 | 600 of 600 |
<!-- /table:evalE_land_basin_e0 -->

<!-- table:evalE_land_basin_full -->
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays monotone to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays monotone to 5 deg, eps 0.05 |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.312 | 3.693 | 1.729 | 45 of 120 | 5.000 | 120 of 120 |
| strict cost, lenient rule | 0.333 | 2.783 | 1.971 | 15 of 120 | 5.000 | 120 of 120 |
| S: sampled binary (control) | 0.348 | 1.688 | 1.315 | 0 of 120 | 5.000 | 111 of 120 |
| E1: distance kernel 3 px | 0.667 | 1.697 | 1.287 | 0 of 120 | 5.000 | 95 of 120 |
| E1: distance kernel 12 px | 1.269 | 1.081 | 0.758 | 0 of 120 | 4.701 | 54 of 120 |
| E2: omega blending, 1 frame | 0.482 | 2.111 | 1.626 | 4 of 120 | 5.000 | 110 of 120 |
| E3: soft peak weights | 0.348 | 1.698 | 1.267 | 0 of 120 | 5.000 | 105 of 120 |
| E4: ray residual, tau 0.05 deg | 0.428 | 5.000 | 5.000 | 117 of 120 | 5.000 | 120 of 120 |
| E4: ray residual, tau 0.25 deg | 1.075 | 5.000 | 3.427 | 87 of 120 | 5.000 | 120 of 120 |
| E1 + E2 | 1.064 | 3.750 | 1.901 | 25 of 120 | 5.000 | 111 of 120 |
| E1 + E2 + E3 | 1.069 | 4.238 | 3.420 | 31 of 120 | 5.000 | 109 of 120 |
| E2 + E3 | 0.481 | 2.228 | 1.787 | 0 of 120 | 5.000 | 103 of 120 |
| E3 + E4 (0.25 deg) | 0.549 | 4.983 | 3.837 | 60 of 120 | 5.000 | 120 of 120 |
<!-- /table:evalE_land_basin_full -->

<!-- table:evalE_land_minima_e0 -->
| cost | rays with a minimum of depth 0.01, dense, 0 to 0.5 deg | rays with a minimum of depth 0.05, dense, 0 to 0.5 deg | rays where the truth is not the lowest point | lowest point when not the truth: median deg (max) | ... and its depth below the truth (median) |
|---|---|---|---|---|---|
| strict cost, C++ rule | 12 of 600 (1.1-3.5%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0005 (0.0010) | 0.0021 |
| strict cost, lenient rule | 13 of 600 (1.3-3.7%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0005 (0.0010) | 0.0021 |
| S: sampled binary (control) | 12 of 600 (1.1-3.5%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E1: distance kernel 3 px | 54 of 600 (7.0-11.6%) | 0 of 600 (0.0-0.6%) | 10 of 600 | 0.0005 (0.0020) | 0.0000 |
| E1: distance kernel 12 px | 61 of 600 (8.0-12.8%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E2: omega blending, 1 frame | 17 of 600 (1.8-4.5%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E3: soft peak weights | 0 of 600 (0.0-0.6%) | 0 of 600 (0.0-0.6%) | 6 of 600 | 0.0008 (0.0030) | 0.0000 |
| E4: ray residual, tau 0.05 deg | 54 of 600 (7.0-11.6%) | 0 of 600 (0.0-0.6%) | 291 of 600 | 0.0075 (0.0625) | 0.0003 |
| E4: ray residual, tau 0.25 deg | 71 of 600 (9.5-14.7%) | 0 of 600 (0.0-0.6%) | 287 of 600 | 0.0080 (0.1535) | 0.0000 |
| E1 + E2 | 88 of 600 (12.1-17.7%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0690 (0.1920) | 0.0080 |
| E1 + E2 + E3 | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 30 of 600 | 0.0075 (0.0235) | 0.0000 |
| E2 + E3 | 0 of 600 (0.0-0.6%) | 0 of 600 (0.0-0.6%) | 6 of 600 | 0.0008 (0.0030) | 0.0000 |
| E3 + E4 (0.25 deg) | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 124 of 600 | 0.0010 (0.0060) | 0.0000 |
<!-- /table:evalE_land_minima_e0 -->

<!-- table:evalE_land_minima_full -->
| cost | rays with a minimum of depth 0.01, dense, 0 to 5 deg | rays with a minimum of depth 0.05, dense, 0 to 5 deg | rays with a minimum of depth 0.01, 0 to 0.5 deg | rays with a minimum of depth 0.05, 0 to 0.5 deg | rays where the truth is not the lowest point | lowest point when not the truth: median deg (max) | ... and its depth below the truth (median) |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 72 of 120 (51.1-68.3%) | 0 of 120 (0.0-3.1%) | 6 of 120 (2.3-10.5%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| strict cost, lenient rule | 97 of 120 (72.9-86.9%) | 0 of 120 (0.0-3.1%) | 5 of 120 (1.8-9.4%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| S: sampled binary (control) | 120 of 120 (96.9-100.0%) | 3 of 120 (0.9-7.1%) | 9 of 120 (4.0-13.6%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E1: distance kernel 3 px | 120 of 120 (96.9-100.0%) | 12 of 120 (5.8-16.7%) | 11 of 120 (5.2-15.7%) | 0 of 120 (0.0-3.1%) | 12 of 120 | 0.0150 (0.0215) | 0.0000 |
| E1: distance kernel 12 px | 120 of 120 (96.9-100.0%) | 43 of 120 (27.8-44.7%) | 9 of 120 (4.0-13.6%) | 0 of 120 (0.0-3.1%) | 37 of 120 | 0.0205 (0.0725) | 0.0000 |
| E2: omega blending, 1 frame | 114 of 120 (89.5-97.7%) | 6 of 120 (2.3-10.5%) | 17 of 120 (9.0-21.5%) | 0 of 120 (0.0-3.1%) | 1 of 120 | 0.2230 (0.2230) | 0.0062 |
| E3: soft peak weights | 120 of 120 (96.9-100.0%) | 7 of 120 (2.9-11.6%) | 5 of 120 (1.8-9.4%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E4: ray residual, tau 0.05 deg | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 54 of 120 | 0.0095 (0.0825) | 0.0002 |
| E4: ray residual, tau 0.25 deg | 31 of 120 (18.8-34.3%) | 0 of 120 (0.0-3.1%) | 24 of 120 (13.8-28.0%) | 0 of 120 (0.0-3.1%) | 55 of 120 | 0.0195 (0.2230) | 0.0001 |
| E1 + E2 | 69 of 120 (48.6-66.0%) | 2 of 120 (0.5-5.9%) | 19 of 120 (10.4-23.4%) | 0 of 120 (0.0-3.1%) | 17 of 120 | 0.1280 (0.2500) | 0.0041 |
| E1 + E2 + E3 | 52 of 120 (34.8-52.3%) | 1 of 120 (0.1-4.6%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E2 + E3 | 119 of 120 (95.4-99.9%) | 10 of 120 (4.6-14.7%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E3 + E4 (0.25 deg) | 52 of 120 (34.8-52.3%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 60 of 120 | 0.0115 (0.0460) | 0.0000 |
<!-- /table:evalE_land_minima_full -->

<!-- table:evalE_land_steps_e0 -->
| cost | cost at the truth (median) | steps per deg, under 0.05 deg | steps per deg, 0.05 to 0.5 deg | steps per deg, 0.5 to 5 deg (4 windows of 0.05 deg) | variation in steps, under 0.05 deg | variation in steps, 0.05 to 0.5 deg | variation in steps, 0.5 to 5 deg (4 windows of 0.05 deg) |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.0606 | 330.8 | 284.7 | 25.8 | 40.7% | 61.9% | 86.2% |
| strict cost, lenient rule | 0.0606 | 318.5 | 349.5 | 50.2 | 38.0% | 55.5% | 84.2% |
| S: sampled binary (control) | 0.0616 | 247.0 | 228.7 | 25.9 | 31.6% | 48.3% | 91.6% |
| E1: distance kernel 3 px | 0.0606 | 42.7 | 42.3 | 13.3 | 94.9% | 64.3% | 52.9% |
| E1: distance kernel 12 px | 0.0606 | 42.7 | 42.3 | 14.2 | 96.0% | 91.4% | 82.8% |
| E2: omega blending, 1 frame | 0.0616 | 217.7 | 234.7 | 36.9 | 20.9% | 37.3% | 68.6% |
| E3: soft peak weights | 0.0000 | 263.2 | 235.0 | 25.3 | 33.4% | 50.3% | 90.9% |
| E4: ray residual, tau 0.05 deg | 0.0977 | 2.0 | 1.4 | 4.9 | 4.9% | 1.2% | 1.3% |
| E4: ray residual, tau 0.25 deg | 0.0049 | 1.8 | 4.3 | 11.1 | 25.6% | 4.8% | 6.4% |
| E1 + E2 | 0.0606 | 2.2 | 1.7 | 1.4 | 46.7% | 3.4% | 1.7% |
| E1 + E2 + E3 | 0.0000 | 0.0 | 0.0 | 0.0 | 0.0% | 0.1% | 0.0% |
| E2 + E3 | 0.0000 | 236.3 | 239.9 | 36.0 | 23.2% | 39.1% | 69.0% |
| E3 + E4 (0.25 deg) | 0.0021 | 249.2 | 229.2 | 34.2 | 31.0% | 40.6% | 27.5% |
<!-- /table:evalE_land_steps_e0 -->

<!-- table:evalE_land_steps_full -->
| cost | cost at the truth (median) | steps per deg, under 0.05 deg | steps per deg, 0.05 to 0.5 deg | steps per deg, 0.5 to 5 deg | variation in steps, under 0.05 deg | variation in steps, 0.05 to 0.5 deg | variation in steps, 0.5 to 5 deg |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.0703 | 113.8 | 261.5 | 25.0 | 88.2% | 84.2% | 94.0% |
| strict cost, lenient rule | 0.0703 | 116.3 | 323.6 | 108.2 | 87.5% | 78.4% | 89.6% |
| S: sampled binary (control) | 0.0717 | 96.8 | 239.0 | 46.9 | 87.6% | 75.1% | 93.5% |
| E1: distance kernel 3 px | 0.0717 | 49.3 | 45.1 | 30.3 | 97.7% | 80.5% | 73.5% |
| E1: distance kernel 12 px | 0.0694 | 49.3 | 49.5 | 47.7 | 97.5% | 93.6% | 89.4% |
| E2: omega blending, 1 frame | 0.0717 | 55.2 | 236.3 | 89.8 | 54.9% | 52.6% | 65.4% |
| E3: soft peak weights | 0.0000 | 96.5 | 241.2 | 45.4 | 87.3% | 75.0% | 92.8% |
| E4: ray residual, tau 0.05 deg | 0.4136 | 0.2 | 1.4 | 8.0 | 0.0% | 1.7% | 2.6% |
| E4: ray residual, tau 0.25 deg | 0.2930 | 4.5 | 4.3 | 13.2 | 0.3% | 6.5% | 9.7% |
| E1 + E2 | 0.0717 | 0.0 | 2.2 | 2.0 | 0.0% | 7.8% | 4.5% |
| E1 + E2 + E3 | 0.0000 | 0.0 | 0.0 | 0.0 | 0.0% | 0.0% | 0.1% |
| E2 + E3 | 0.0000 | 55.8 | 237.9 | 88.1 | 56.3% | 52.9% | 65.5% |
| E3 + E4 (0.25 deg) | 0.1464 | 96.3 | 239.7 | 54.8 | 84.2% | 63.0% | 55.5% |
<!-- /table:evalE_land_steps_full -->

<!-- table:evalE_land_types_e0 -->
| cost | isolated steps | at a frame event | at a pixel event | at a presence event | at a detector event | at no strict event |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 89902 | 10.4% | 89.1% | 0.4% | 0.1% | 0.0% |
| strict cost, lenient rule | 109934 | 10.0% | 89.6% | 0.3% | 0.1% | 0.0% |
| S: sampled binary (control) | 72265 | 13.4% | 82.0% | 0.5% | 0.2% | 4.0% |
| E1: distance kernel 3 px | 14299 | 93.7% | 1.5% | 3.4% | 1.1% | 0.2% |
| E1: distance kernel 12 px | 14414 | 93.6% | 1.6% | 3.5% | 1.2% | 0.2% |
| E2: omega blending, 1 frame | 74322 | 4.1% | 90.5% | 0.5% | 0.2% | 4.7% |
| E3: soft peak weights | 74390 | 13.2% | 82.7% | 0.1% | 0.1% | 3.9% |
| E4: ray residual, tau 0.05 deg | 1014 | 41.6% | 14.2% | 42.8% | 0.3% | 1.1% |
| E4: ray residual, tau 0.25 deg | 2540 | 74.3% | 6.8% | 18.1% | 0.1% | 0.6% |
| E1 + E2 | 696 | 2.2% | 3.2% | 68.0% | 26.7% | 0.0% |
| E1 + E2 + E3 | 8 | 0.0% | 12.5% | 87.5% | 0.0% | 0.0% |
| E2 + E3 | 76196 | 4.2% | 91.0% | 0.1% | 0.1% | 4.6% |
| E3 + E4 (0.25 deg) | 73457 | 15.1% | 80.8% | 0.1% | 0.1% | 3.9% |
<!-- /table:evalE_land_types_e0 -->

<!-- table:evalE_land_types_full -->
| cost | isolated steps | at a frame event | at a pixel event | at a presence event | at a detector event | at no strict event |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 28290 | 15.2% | 83.1% | 1.6% | 0.1% | 0.0% |
| strict cost, lenient rule | 76615 | 10.6% | 88.5% | 0.8% | 0.1% | 0.0% |
| S: sampled binary (control) | 38789 | 18.5% | 76.5% | 1.7% | 0.7% | 2.6% |
| E1: distance kernel 3 px | 19069 | 92.2% | 2.1% | 3.6% | 1.8% | 0.2% |
| E1: distance kernel 12 px | 28709 | 92.9% | 2.2% | 2.6% | 2.1% | 0.2% |
| E2: omega blending, 1 frame | 61571 | 4.0% | 91.3% | 1.0% | 0.5% | 3.2% |
| E3: soft peak weights | 38094 | 18.8% | 78.3% | 0.2% | 0.1% | 2.6% |
| E4: ray residual, tau 0.05 deg | 4409 | 78.1% | 6.8% | 14.0% | 0.4% | 0.7% |
| E4: ray residual, tau 0.25 deg | 7400 | 86.4% | 4.6% | 8.3% | 0.2% | 0.5% |
| E1 + E2 | 1189 | 2.0% | 9.2% | 54.5% | 33.6% | 0.8% |
| E1 + E2 + E3 | 19 | 5.3% | 42.1% | 52.6% | 0.0% | 0.0% |
| E2 + E3 | 60777 | 4.0% | 92.6% | 0.1% | 0.1% | 3.2% |
| E3 + E4 (0.25 deg) | 43125 | 28.8% | 68.6% | 0.2% | 0.1% | 2.3% |
<!-- /table:evalE_land_types_full -->

<!-- table:evalE_rank_e0_cpp_1deg -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 109 of 111 | 0 | 142 of 147 | 0 | 144 of 147 | 0 | 395 of 405 (95.5-98.7%) |  |
| S | 109 of 111 (+0 / -0) | 0 | 144 of 147 (+2 / -0) | 0 | 144 of 147 (+0 / -0) | 0 | 397 of 405 (96.2-99.0%) | 2 / 0 |
| E1s3 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E1s8 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E1a | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E2a | 110 of 111 (+1 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 143 of 147 (+0 / -1) | 0 | 400 of 405 (97.1-99.5%) | 6 / 1 |
| E3 | 111 of 111 (+2 / -0) | 0 | 142 of 147 (+2 / -2) | 0 | 144 of 147 (+0 / -0) | 0 | 397 of 405 (96.2-99.0%) | 4 / 2 |
| E4t05 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 146 of 147 (+2 / -0) | 0 | 404 of 405 (98.6-100.0%) | 9 / 0 |
| E4t25 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E4a | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E12a | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E123a | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E23a | 111 of 111 (+2 / -0) | 0 | 146 of 147 (+5 / -1) | 0 | 144 of 147 (+1 / -1) | 0 | 401 of 405 (97.5-99.6%) | 8 / 2 |
| E34a | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E1234a | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E1s12 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E2w1 | 110 of 111 (+1 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 143 of 147 (+0 / -1) | 0 | 400 of 405 (97.1-99.5%) | 6 / 1 |
| E12s3 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E123s3 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
| E23w1 | 111 of 111 (+2 / -0) | 0 | 146 of 147 (+5 / -1) | 0 | 144 of 147 (+1 / -1) | 0 | 401 of 405 (97.5-99.6%) | 8 / 2 |
| E34t25 | 111 of 111 (+2 / -0) | 0 | 147 of 147 (+5 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 405 of 405 (99.1-100.0%) | 10 / 0 |
<!-- /table:evalE_rank_e0_cpp_1deg -->

<!-- table:evalE_rank_e0_cpp_next -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 181 of 196 | 3 | 149 of 165 | 0 | 144 of 147 | 0 | 474 of 508 (90.8-95.2%) |  |
| S | 184 of 196 (+5 / -2) | 2 | 150 of 165 (+2 / -1) | 0 | 144 of 147 (+0 / -0) | 0 | 478 of 508 (91.7-95.8%) | 7 / 3 |
| E1s3 | 195 of 196 (+14 / -0) | 0 | 163 of 165 (+14 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 505 of 508 (98.3-99.8%) | 31 / 0 |
| E1s8 | 195 of 196 (+14 / -0) | 0 | 163 of 165 (+14 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 505 of 508 (98.3-99.8%) | 31 / 0 |
| E1a | 195 of 196 (+14 / -0) | 0 | 163 of 165 (+14 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 505 of 508 (98.3-99.8%) | 31 / 0 |
| E2a | 186 of 196 (+8 / -3) | 1 | 159 of 165 (+10 / -0) | 0 | 143 of 147 (+0 / -1) | 0 | 488 of 508 (94.0-97.4%) | 18 / 4 |
| E3 | 177 of 196 (+4 / -8) | 2 | 149 of 165 (+4 / -4) | 0 | 144 of 147 (+0 / -0) | 0 | 470 of 508 (89.9-94.5%) | 8 / 12 |
| E4t05 | 193 of 196 (+12 / -0) | 0 | 163 of 165 (+14 / -0) | 0 | 146 of 147 (+2 / -0) | 0 | 502 of 508 (97.4-99.5%) | 28 / 0 |
| E4t25 | 195 of 196 (+14 / -0) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 507 of 508 (98.9-100.0%) | 33 / 0 |
| E4a | 196 of 196 (+15 / -0) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 508 of 508 (99.2-100.0%) | 34 / 0 |
| E12a | 196 of 196 (+15 / -0) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 508 of 508 (99.2-100.0%) | 34 / 0 |
| E123a | 196 of 196 (+15 / -0) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 508 of 508 (99.2-100.0%) | 34 / 0 |
| E23a | 186 of 196 (+8 / -3) | 1 | 156 of 165 (+9 / -2) | 0 | 144 of 147 (+1 / -1) | 0 | 486 of 508 (93.5-97.1%) | 18 / 6 |
| E34a | 195 of 196 (+14 / -0) | 0 | 164 of 165 (+15 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 506 of 508 (98.6-99.9%) | 32 / 0 |
| E1234a | 196 of 196 (+15 / -0) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 508 of 508 (99.2-100.0%) | 34 / 0 |
| E1s12 | 195 of 196 (+14 / -0) | 0 | 163 of 165 (+14 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 505 of 508 (98.3-99.8%) | 31 / 0 |
| E2w1 | 184 of 196 (+8 / -5) | 1 | 159 of 165 (+10 / -0) | 0 | 143 of 147 (+0 / -1) | 0 | 486 of 508 (93.5-97.1%) | 18 / 6 |
| E12s3 | 195 of 196 (+14 / -0) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 507 of 508 (98.9-100.0%) | 33 / 0 |
| E123s3 | 194 of 196 (+14 / -1) | 0 | 165 of 165 (+16 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 506 of 508 (98.6-99.9%) | 33 / 1 |
| E23w1 | 185 of 196 (+8 / -4) | 1 | 156 of 165 (+9 / -2) | 0 | 144 of 147 (+1 / -1) | 0 | 485 of 508 (93.3-97.0%) | 18 / 7 |
| E34t25 | 194 of 196 (+13 / -0) | 0 | 164 of 165 (+15 / -0) | 0 | 147 of 147 (+3 / -0) | 0 | 505 of 508 (98.3-99.8%) | 31 / 0 |
<!-- /table:evalE_rank_e0_cpp_next -->

<!-- table:evalE_rank_e0_old_1deg -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 115 of 115 | 0 | 158 of 160 | 0 | 159 of 161 | 0 | 432 of 436 (97.7-99.6%) |  |
| S | 115 of 115 (+0 / -0) | 0 | 156 of 160 (+0 / -2) | 0 | 159 of 161 (+0 / -0) | 0 | 430 of 436 (97.0-99.4%) | 0 / 2 |
| E1s3 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E1s8 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E1a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E2a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 160 of 161 (+1 / -0) | 0 | 435 of 436 (98.7-100.0%) | 3 / 0 |
| E3 | 115 of 115 (+0 / -0) | 0 | 157 of 160 (+0 / -1) | 0 | 159 of 161 (+0 / -0) | 0 | 431 of 436 (97.3-99.5%) | 0 / 1 |
| E4t05 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E4t25 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E4a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E12a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E123a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E23a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 160 of 161 (+1 / -0) | 0 | 435 of 436 (98.7-100.0%) | 3 / 0 |
| E34a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E1234a | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E1s12 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E2w1 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 160 of 161 (+1 / -0) | 0 | 435 of 436 (98.7-100.0%) | 3 / 0 |
| E12s3 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E123s3 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
| E23w1 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 160 of 161 (+1 / -0) | 0 | 435 of 436 (98.7-100.0%) | 3 / 0 |
| E34t25 | 115 of 115 (+0 / -0) | 0 | 160 of 160 (+2 / -0) | 0 | 161 of 161 (+2 / -0) | 0 | 436 of 436 (99.1-100.0%) | 4 / 0 |
<!-- /table:evalE_rank_e0_old_1deg -->

<!-- table:evalE_rank_e0_old_next -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 186 of 197 | 1 | 166 of 176 | 0 | 160 of 163 | 0 | 512 of 536 (93.4-97.0%) |  |
| S | 188 of 197 (+3 / -1) | 1 | 164 of 176 (+0 / -2) | 0 | 160 of 163 (+0 / -0) | 0 | 512 of 536 (93.4-97.0%) | 3 / 3 |
| E1s3 | 195 of 197 (+9 / -0) | 0 | 175 of 176 (+9 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 532 of 536 (98.1-99.7%) | 20 / 0 |
| E1s8 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E1a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E2a | 188 of 197 (+4 / -2) | 0 | 169 of 176 (+3 / -0) | 0 | 160 of 163 (+1 / -1) | 0 | 517 of 536 (94.5-97.7%) | 8 / 3 |
| E3 | 188 of 197 (+4 / -2) | 1 | 163 of 176 (+1 / -4) | 0 | 159 of 163 (+0 / -1) | 0 | 510 of 536 (93.0-96.7%) | 5 / 7 |
| E4t05 | 195 of 197 (+10 / -1) | 0 | 174 of 176 (+8 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 531 of 536 (97.8-99.6%) | 20 / 1 |
| E4t25 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E4a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E12a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E123a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E23a | 188 of 197 (+4 / -2) | 0 | 170 of 176 (+5 / -1) | 0 | 160 of 163 (+1 / -1) | 0 | 518 of 536 (94.8-97.9%) | 10 / 4 |
| E34a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E1234a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E1s12 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E2w1 | 186 of 197 (+2 / -2) | 1 | 167 of 176 (+3 / -2) | 0 | 160 of 163 (+1 / -1) | 0 | 513 of 536 (93.6-97.1%) | 6 / 5 |
| E12s3 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E123s3 | 196 of 197 (+10 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E23w1 | 186 of 197 (+3 / -3) | 1 | 169 of 176 (+5 / -2) | 0 | 160 of 163 (+1 / -1) | 0 | 515 of 536 (94.1-97.4%) | 9 / 6 |
| E34t25 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
<!-- /table:evalE_rank_e0_old_next -->

<!-- table:evalE_rank_full_cpp_1deg -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 3 of 3 | 0 | 6 of 6 | 0 | 6 of 6 | 0 | 15 of 15 (79.6-100.0%) |  |
| S | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E1s3 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E1s8 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E1a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E2a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E3 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E4t05 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E4t25 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E4a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E12a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E123a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E23a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E34a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E1234a | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E1s12 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E2w1 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E12s3 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E123s3 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E23w1 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
| E34t25 | 3 of 3 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 15 of 15 (79.6-100.0%) | 0 / 0 |
<!-- /table:evalE_rank_full_cpp_1deg -->

<!-- table:evalE_rank_full_cpp_next -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 6 of 6 | 2 | 6 of 6 | 0 | 6 of 6 | 0 | 18 of 18 (82.4-100.0%) |  |
| S | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1s3 | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1s8 | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1a | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E2a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E3 | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E4t05 | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E4t25 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E4a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E12a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E123a | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E23a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E34a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1234a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1s12 | 6 of 6 (+0 / -0) | 1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E2w1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E12s3 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E123s3 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E23w1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E34t25 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
<!-- /table:evalE_rank_full_cpp_next -->

<!-- table:evalE_rank_full_old_1deg -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 4 of 4 | 0 | 6 of 6 | 0 | 6 of 6 | 0 | 16 of 16 (80.6-100.0%) |  |
| S | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E1s3 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E1s8 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E1a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E2a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E3 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E4t05 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E4t25 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E4a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E12a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E123a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E23a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E34a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E1234a | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E1s12 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E2w1 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E12s3 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E123s3 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E23w1 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
| E34t25 | 4 of 4 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 16 of 16 (80.6-100.0%) | 0 / 0 |
<!-- /table:evalE_rank_full_old_1deg -->

<!-- table:evalE_rank_full_old_next -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 6 of 6 | 1 | 6 of 6 | 0 | 6 of 6 | 0 | 18 of 18 (82.4-100.0%) |  |
| S | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1s3 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1s8 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E2a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E3 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E4t05 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E4t25 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E4a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E12a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E123a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E23a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E34a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1234a | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E1s12 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E2w1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E12s3 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E123s3 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E23w1 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
| E34t25 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 6 of 6 (+0 / -0) | 0 | 18 of 18 (82.4-100.0%) | 0 / 0 |
<!-- /table:evalE_rank_full_old_next -->

<!-- table:evalE_timing_ms -->
| cost | full-sample: ms per evaluation, batch of 256 | full-sample: ms per evaluation, one per call | full-sample: x strict | E0: ms per evaluation, batch of 256 | E0: ms per evaluation, one per call | E0: x strict |
|---|---|---|---|---|---|---|
| strict cost (one call per orientation) | 0.381 | 0.381 | 1.000 | 0.392 | 0.392 | 1.000 |
| S | 0.470 | 0.809 | 1.233 | 0.228 | 0.601 | 0.581 |
| E1 | 1.332 | 1.649 | 3.492 | 0.240 | 0.634 | 0.613 |
| E2 | 2.080 | 2.675 | 5.453 | 0.886 | 1.590 | 2.261 |
| E3 | 0.473 | 0.844 | 1.240 | 0.227 | 0.640 | 0.580 |
| E4 | 1.186 | 2.153 | 3.109 | 0.313 | 1.248 | 0.798 |
| E1+E2 | 6.663 | 7.268 | 17.466 | 0.950 | 1.710 | 2.424 |
| E1+E2+E3 | 6.728 | 7.485 | 17.637 | 0.952 | 1.747 | 2.430 |
| E2+E3 | 2.136 | 2.792 | 5.599 | 0.888 | 1.627 | 2.265 |
| E3+E4 | 1.639 | 2.811 | 4.297 | 0.495 | 1.632 | 1.262 |
| E1+E2+E3+E4 | 7.903 | 9.499 | 20.718 | 1.223 | 2.745 | 3.120 |
<!-- /table:evalE_timing_ms -->

<!-- table:evalE_timing_per_voxel -->
| voxel set | arm | strict cost time per voxel, baseline (s) | strict cost time per voxel, arm (s) | smooth evaluations per voxel | smooth cost time per voxel (s) | one-off build per image set (s) | extra time with a build per voxel | extra time with one build per image set |
|---|---|---|---|---|---|---|---|---|
| E0 200 | E3+E4,tau=0.25 | 19.1 | 19.3 | 252 | 0.1 | 3.6 | +20.7% | +1.7% |
| six full-sample | E3+E4,tau=0.25 | 126.7 | 127.6 | 15030 | 24.6 | 16.6 | +33.3% | +20.2% |
| E0 200 | E1+E2+E3+E4 | 19.1 | 19.5 | 263 | 0.3 | 3.6 | +22.6% | +3.6% |
| six full-sample | E1+E2+E3+E4 | 126.7 | 132.8 | 15755 | 124.5 | 16.6 | +116.2% | +103.1% |
<!-- /table:evalE_timing_per_voxel -->
