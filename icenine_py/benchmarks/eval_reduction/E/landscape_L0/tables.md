<!-- table:evalE_landL0_basin_e0 -->
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays with no drop over eps out to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays with no drop over eps out to 5 deg, eps 0.05 |
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
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays with no drop over eps out to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays with no drop over eps out to 5 deg, eps 0.05 |
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
