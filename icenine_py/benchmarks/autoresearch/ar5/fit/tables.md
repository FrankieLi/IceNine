<!-- table:base -->
| item | new | stored |
|---|---|---|
| evaluations before the pass | 7201863 | 7201863 |
| wrong voxels before the pass | 67 | 67 |
| seed steps | 329 | 329 |
<!-- /table:base -->

<!-- table:err -->
| sample | arm | steps with candidates | median error (deg) | 90th percentile (deg) | within 0.2 deg (%) |
|---|---|---|---|---|---|
| big | M | 329 | 0.016 | 0.056 | 100.0 |
| big | C60 | 329 | 0.126 | 0.542 | 68.7 |
| big | C200 | 329 | 0.013 | 0.050 | 97.0 |
| big | C0 | 329 | 0.453 | 38.956 | 13.1 |
| first gradient sample | M | 1279 | 0.017 | 0.080 | 95.9 |
| first gradient sample | C60 | 1279 | 0.058 | 17.465 | 74.2 |
| first gradient sample | C200 | 1279 | 0.010 | 0.234 | 87.3 |
| first gradient sample | C0 | 1279 | 0.276 | 44.539 | 40.3 |
| mid | M | 1611 | 0.020 | 0.076 | 97.8 |
| mid | C60 | 1611 | 0.160 | 23.341 | 56.0 |
| mid | C200 | 1611 | 0.013 | 0.317 | 85.5 |
| mid | C0 | 1611 | 0.423 | 41.431 | 16.1 |
<!-- /table:err -->

<!-- table:fits -->
| sample | S8 fits | S8 MC fallbacks | S8 MC ended lower | S10 fits | S10 CMA cost did not fall | S8 neighbour + revisit evaluations | S10 neighbour + revisit evaluations | S8 seed steps | S10 seed steps |
|---|---|---|---|---|---|---|---|---|---|
| big | 32563 | 3197 | 3 | 32398 | 3152 | 4345151 | 1976278 | 329 | 333 |
| first gradient sample | 53300 | 3903 | 3 | 53386 | 3876 | 5858731 | 3256546 | 1279 | 1327 |
| mid | 51696 | 4568 | 45 | 52082 | 4496 | 6376933 | 3177002 | 1611 | 1604 |
<!-- /table:fits -->

<!-- table:mid_pass -->
| arm | rule (a) | rule (b) | evaluations fall (%) | rule (c) | passes |
|---|---|---|---|---|---|
| C60 | False | False | 96.3 | True | False |
| C200 | False | False | 93.6 | True | False |
| C0 | False | False | 99.7 | True | False |
<!-- /table:mid_pass -->

<!-- table:part0 -->
| sample | seed steps | refinements | refinement evaluations | mean per refinement | MC fallbacks | share of fits (%) | fallback evaluations | mean per fallback | ended lower than CMA | Wilson 95% interval | evaluations of S8 with the pass | refinements + fallbacks (% of total) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| big | 329 | 465 | 2850132 | 6129 | 3197 | 9.8 | 2358808 | 738 | 3 | 0.03 to 0.28% | 8340428 | 62.5 |
| first gradient sample | 1279 | 2356 | 2983177 | 1266 | 3903 | 7.3 | 2607431 | 668 | 3 | 0.03 to 0.23% | 18000230 | 31.1 |
| mid | 1611 | 3217 | 11975132 | 3722 | 4568 | 8.8 | 3223477 | 706 | 45 | 0.74 to 1.32% | 26443074 | 57.5 |
<!-- /table:part0 -->

<!-- table:part1 -->
| sample | arm | steps | pass the seed test | Wilson 95% interval | passes: M only | passes: arm only | within 1 degree | within 1 deg: M only | within 1 deg: arm only | seed-step evaluations | fewer than M (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| big | M | 329 | 329 | 98.8 to 100.0% | 0 | 0 | 329 | 0 | 0 | 2856712 | 0.0 |
| big | C60 | 329 | 244 | 69.2 to 78.6% | 85 | 0 | 303 | 26 | 0 | 88503 | 96.9 |
| big | C200 | 329 | 323 | 96.1 to 99.2% | 6 | 0 | 324 | 5 | 0 | 117733 | 95.9 |
| big | C0 | 329 | 43 | 9.8 to 17.1% | 286 | 0 | 254 | 75 | 0 | 6580 | 99.8 |
| first gradient sample | M | 1279 | 1263 | 98.0 to 99.2% | 0 | 0 | 1254 | 0 | 0 | 3008757 | 0.0 |
| first gradient sample | C60 | 1279 | 1103 | 84.2 to 88.0% | 162 | 2 | 1145 | 111 | 2 | 273789 | 90.9 |
| first gradient sample | C200 | 1279 | 1209 | 93.1 to 95.6% | 58 | 4 | 1222 | 40 | 8 | 532703 | 82.3 |
| first gradient sample | C0 | 1279 | 499 | 36.4 to 41.7% | 764 | 0 | 923 | 335 | 4 | 25580 | 99.1 |
| mid | M | 1611 | 1596 | 98.5 to 99.4% | 0 | 0 | 1595 | 0 | 0 | 12007352 | 0.0 |
| mid | C60 | 1611 | 1037 | 62.0 to 66.7% | 559 | 0 | 1445 | 154 | 4 | 447142 | 96.3 |
| mid | C200 | 1611 | 1467 | 89.6 to 92.4% | 133 | 4 | 1542 | 58 | 5 | 765267 | 93.6 |
| mid | C0 | 1611 | 180 | 9.7 to 12.8% | 1416 | 0 | 1280 | 319 | 4 | 32220 | 99.7 |
<!-- /table:part1 -->

<!-- table:part2 -->
| run | sample | seed step | grid fallback | neighbour fits | revisits | boundary pass | total | grid fallbacks (count) | wrong before the pass | wrong after the pass | fixed / broken before | fixed / broken after | median error after, all 24570 voxels (deg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S8 | big | 2856712 | 0 | 3491437 | 853714 | 1138565 | 8340428 | 0 | 67 (0.21 to 0.35%) | 24 (0.07 to 0.15%) |  |  | 0.0316 |
| S8 | first gradient sample | 3008757 | 5449501 | 2798001 | 3060730 | 3683241 | 18000230 | 16 | 210 (0.75 to 0.98%) | 63 (0.20 to 0.33%) |  |  | 0.0226 |
| S8 | mid | 12007352 | 5002750 | 3859661 | 2517272 | 3056039 | 26443074 | 15 | 84 (0.28 to 0.42%) | 31 (0.09 to 0.18%) |  |  | 0.0297 |
| S10 | big | 2790857 | 0 | 1496574 | 479704 | 1131123 | 5898258 | 0 | 77 (0.25 to 0.39%) | 22 (0.06 to 0.14%) | 10 / 20 | 11 / 9 | 0.0318 |
| S10 | first gradient sample | 3122117 | 5108107 | 1477481 | 1779065 | 3681228 | 15167998 | 15 | 198 (0.70 to 0.93%) | 67 (0.21 to 0.35%) | 51 / 39 | 29 / 33 | 0.0224 |
| S10 | mid | 11431628 | 5008366 | 1489254 | 1687748 | 3051708 | 22668704 | 15 | 73 (0.24 to 0.37%) | 35 (0.10 to 0.20%) | 25 / 14 | 11 / 15 | 0.0300 |
<!-- /table:part2 -->

<!-- table:pick -->
| arm | pass the seed test | within 1 degree | seed-step evaluations | rule (a) | rule (b) | picked |
|---|---|---|---|---|---|---|
| M | 1592 | 1583 | 5865469 | reference | reference | yes |
| C60 | 1347 | 1448 | 362292 | False | False | no |
| C200 | 1532 | 1546 | 650436 | False | False | no |
| C0 | 542 | 1177 | 32160 | False | False | no |
<!-- /table:pick -->

<!-- table:pred -->
| prediction | observed | verdict |
|---|---|---|
| 1 | big 3 of 3197 (0.09%); first gradient sample 3 of 3903 (0.08%); mid 45 of 4568 (0.99%) | agrees |
| 2 | pick M; C60 cuts the seed step by big 96.9%, first gradient sample 90.9%, mid 96.3% | no (the pick is M); the cut holds |
| 3 | C0 passes the seed test in 542 of 1608 steps, M in 1592 | agrees |
| 4 | big 5898258 (limit 4000000; rule passes True); first gradient sample 15167998 (limit 14000000; rule passes False); mid 22668704 (limit 13500000; rule passes False) | no |
| 5 | big 0 against 0; first gradient sample 15 against 16; mid 15 against 15 | agrees |
<!-- /table:pred -->

<!-- table:rule -->
| run | sample | fewer than S8 (%) | rule (a) | broken minus fixed | rule (b) | median error minus S8's, all 24570 voxels (deg) | rule (c) | passes |
|---|---|---|---|---|---|---|---|---|
| S10 | big | 29.3 | True | -2 | True | +0.0002 | True | True |
| S10 | first gradient sample | 15.7 | False | 4 | True | -0.0002 | True | False |
| S10 | mid | 14.3 | False | 4 | True | +0.0002 | True | False |
<!-- /table:rule -->

<!-- table:shares -->
| sample | S8 total | S10 total | saved | S8 seed step (% of total) | S8 grid fallback (%) | S8 neighbour + revisit (%) | S8 boundary pass (%) | saved in neighbour + revisit (% of S8 total) |
|---|---|---|---|---|---|---|---|---|
| big | 8340428 | 5898258 | 2442170 | 34.3 | 0.0 | 52.1 | 13.7 | 28.4 |
| first gradient sample | 18000230 | 15167998 | 2832232 | 16.7 | 30.3 | 32.5 | 20.5 | 14.5 |
| mid | 26443074 | 22668704 | 3774370 | 45.4 | 18.9 | 24.1 | 11.6 | 12.1 |
<!-- /table:shares -->
