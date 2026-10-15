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
