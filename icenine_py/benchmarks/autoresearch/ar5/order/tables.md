<!-- table:order_accuracy -->
| sample | arm | wb | wa | med | f8 | fewer | r5 | wall |
|---|---|---|---|---|---|---|---|---|
| held | A | 104 (0.35 to 0.51%) | 73 (0.24 to 0.37%) | 0.02607 | 25 / 28 | 22.4% | pass | 127 |
| held | Y | 104 (0.35 to 0.51%) | 78 (0.25 to 0.40%) | 0.02603 | 23 / 31 | 28.3% | fail | 123 |
| mid | A | 75 (0.24 to 0.38%) | 35 (0.10 to 0.20%) | 0.02862 | 9 / 13 | 40.8% | pass | 39 |
| mid | Y | 76 (0.25 to 0.39%) | 38 (0.11 to 0.21%) | 0.02861 | 7 / 14 | 48.8% | fail | 38 |
| grad | A | 206 (0.73 to 0.96%) | 56 (0.18 to 0.30%) | 0.02135 | 31 / 24 | 32.3% | pass | 30 |
| grad | Y | 207 (0.74 to 0.96%) | 56 (0.18 to 0.30%) | 0.02135 | 31 / 24 | 33.7% | pass | 30 |
| big | A | 74 (0.24 to 0.38%) | 21 (0.06 to 0.13%) | 0.03044 | 8 / 5 | 54.3% | pass | 9 |
| big | Y | 72 (0.23 to 0.37%) | 18 (0.05 to 0.12%) | 0.03040 | 11 / 5 | 59.2% | pass | 9 |
| big | Yl | 72 (0.23 to 0.37%) | 25 (0.07 to 0.15%) | 0.03071 | 8 / 9 | 59.0% | pass | 9 |
<!-- /table:order_accuracy -->

<!-- table:order_predictions -->
| sample | p12 | p3 | p4 | p5 |
|---|---|---|---|---|
| held | met (seed 22.6%, total 7.7%; range 20-45 / 6-16) | not met (grid 65 to 64, net 5) | met (-0.00004) | met (0 of 10564 shared fits differ) |
| mid | met (seed 38.8%, total 13.4%; range 30-45 / 11-18) | met (grid 12 to 12, net 3) | met (-0.00001) | met (0 of 2896 shared fits differ) |
| grad | met (seed 16.4%, total 2.1%; range 10-25 / 1-5) | met (grid 11 to 11, net 0) | met (+0.00000) | met (0 of 2649 shared fits differ) |
| big | not met (seed 60.1%, total 10.6%; range 60-80 / 12-22) | met (grid 0 to 0, net -3) | met (-0.00004) | met (0 of 318 shared fits differ) |
<!-- /table:order_predictions -->

<!-- table:order_rules -->
| sample | r1 | r2 | r3 | r4 | r5 | a5 | y |
|---|---|---|---|---|---|---|---|
| held | pass | pass (5) | pass | pass | fail (28.3% fewer, net 8, median -0.00256) | pass (22.4% fewer, net 3) | fail |
| mid | pass | pass (3) | pass | pass | fail (48.8% fewer, net 7, median -0.00112) | pass (40.8% fewer, net 4) | fail |
| grad | pass | pass (0) | pass | pass | pass (33.7% fewer, net -7, median -0.00128) | pass (32.3% fewer, net -7) | pass |
| big | pass | pass (-3) | pass | pass | pass (59.2% fewer, net -6, median -0.00121) | pass (54.3% fewer, net -3) | pass |
<!-- /table:order_rules -->

<!-- table:order_stages -->
| sample | arm | steps | seed_evals | mc | mc_steps | nb | rv | grid_e | before | pas | total | grid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| held | A | 3709 | 14,173,739 | 5457 | 2234 | 1,480,653 | 3,211,833 | 22,458,034 | 41,324,259 | 5,192,137 | 46,516,396 | 65 |
| held | Y | 3680 | 10,975,888 | 4100 | 1300 | 1,480,531 | 3,190,788 | 22,112,176 | 37,759,383 | 5,190,246 | 42,949,629 | 64 |
| mid | A | 1607 | 5,417,024 | 1605 | 770 | 1,488,949 | 1,678,659 | 4,007,241 | 12,591,873 | 3,053,599 | 15,645,472 | 12 |
| mid | Y | 1610 | 3,316,028 | 1047 | 342 | 1,488,949 | 1,684,210 | 4,007,241 | 10,496,428 | 3,050,610 | 13,547,038 | 12 |
| grad | A | 1278 | 1,527,842 | 991 | 478 | 1,477,542 | 1,778,394 | 3,726,723 | 8,510,501 | 3,679,276 | 12,189,777 | 11 |
| grad | Y | 1270 | 1,277,801 | 662 | 225 | 1,477,542 | 1,774,673 | 3,726,723 | 8,256,739 | 3,679,337 | 11,936,076 | 11 |
| big | A | 317 | 687,902 | 153 | 99 | 1,496,513 | 488,000 | 0 | 2,672,415 | 1,135,332 | 3,807,747 | 0 |
| big | Y | 319 | 274,195 | 46 | 20 | 1,496,391 | 495,503 | 0 | 2,266,089 | 1,138,138 | 3,404,227 | 0 |
| big | Yl | 345 | 300,088 | 70 | 26 | 1,496,696 | 488,305 | 0 | 2,285,089 | 1,138,260 | 3,423,349 | 0 |
<!-- /table:order_stages -->

<!-- table:order_y_vs_a -->
| sample | seed | sc | tc | fb | fbb | wrong | grid | med | dif |
|---|---|---|---|---|---|---|---|---|---|
| held | 14,173,739 to 10,975,888 | 22.6% | 7.7% | 1 / 6 | 4 / 4 | 73 to 78 | 65 to 64 | 0.02607 to 0.02603 | 232 of 3231 (0 wrong in A, 0 in Y) |
| mid | 5,417,024 to 3,316,028 | 38.8% | 13.4% | 2 / 5 | 1 / 2 | 35 to 38 | 12 to 12 | 0.02862 to 0.02861 | 115 of 1173 (0 wrong in A, 0 in Y) |
| grad | 1,527,842 to 1,277,801 | 16.4% | 2.1% | 0 / 0 | 0 / 1 | 56 to 56 | 11 to 11 | 0.02135 to 0.02135 | 5 of 1259 (0 wrong in A, 0 in Y) |
| big | 687,902 to 274,195 | 60.1% | 10.6% | 3 / 0 | 2 / 0 | 21 to 18 | 0 to 0 | 0.03044 to 0.03040 | 22 of 217 (0 wrong in A, 0 in Y) |
<!-- /table:order_y_vs_a -->

<!-- table:order_yl -->
| what | Y | Yl |
|---|---|---|
| seed-step evaluations | 274,195 | 300,088 |
| total with the pass | 3,404,227 | 3,423,349 |
| wrong after the pass | 18 | 25 |
| grid fallbacks | 0 | 0 |
| fixed / broken, Yl against Y (after the pass) |  | 5 / 12 |
| fixed / broken, Yl against A (after the pass) |  | 5 / 9 |
<!-- /table:order_yl -->
