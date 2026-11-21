<!-- table:broken -->
| broken voxel | error in J5 (degree) | error in J5L (degree) | fits skipped there | fits run there |
|---|---|---|---|---|
| 1091 | 0.010 | 44.9 | 2 | 8 |
| 1426 | 0.075 | 30.3 | 0 | 1 |
| 2078 | 0.022 | 24.2 | 0 | 1 |
| 3324 | 0.068 | 43.7 | 0 | 5 |
<!-- /table:broken -->

<!-- table:judge_edit -->
| rule on mid | frozen grid comparison | edited grid comparison |
|---|---|---|
| (1) 3% fewer than J5 | True | True |
| (2) broken minus fixed against J5 at most 5 | True | True |
| (3) median at most 0.001 degree higher | True | True |
| (4) wrong seeds, grid searches (2 against 2), reseeds | True | True |
| (5) against S8 (J5 passes it) | False | False |
| J5L passes on mid | False | False |
<!-- /table:judge_edit -->

<!-- table:live_mid -->
| arm (mid, with the pass) | evaluations | broken / fixed against J5 | broken / fixed against S8 | median error (degree) |
|---|---|---|---|---|
| J5L, live, gate model | 9,990,663 | 4 / 0 | 17 / 9 | 0.02857 |
| J5 | 10,589,493 |  | 13 / 9 | 0.02860 |
| S8 | 26,443,074 |  |  | 0.02972 |
<!-- /table:live_mid -->

<!-- table:loso -->
| held-out sample | net saving | wrong skips (allowed) | changed final | S(0.3, both) | ceiling |
|---|---|---|---|---|---|
| held3 | 5.46% | 30 (72) | 1 | 3.98%, 21 wrong | 13.31% |
| held | 6.16% | 70 (76) | 0 | 4.21%, 34 wrong | 14.90% |
| mid | 5.79% | 46 (52) | 0 | 3.62%, 26 wrong | 16.75% |
| grad | 11.83% | 43 (53) | 0 | 8.74%, 24 wrong | 20.31% |
| big | 5.38% | 17 (32) | 0 | 3.20%, 6 wrong | 14.19% |
<!-- /table:loso -->

<!-- table:replay_vs_live -->
| quantity (mid) | replay | live |
|---|---|---|
| fits | 52,027 | 52,075 |
| saving as a share of J5 with the pass | 5.79% | 5.65% |
| skipped fits | 10,949 | 10,836 |
| fits compared up to the replay's first wrong skip | 7,232 | 7,232 |
| of these, fits with a feature difference (feature 13 only) | 0 | 240 |
| first fit with a difference | none | 1016 |
<!-- /table:replay_vs_live -->
