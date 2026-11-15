<!-- table:repair_accuracy -->
| sample | wrong_before | wrong_after | fixed | broken | med_before | med_after |
|---|---|---|---|---|---|---|
| held3 | 54 | 51 | 3 | 0 | 0.02545 | 0.02543 |
| held | 67 | 58 | 9 | 0 | 0.02596 | 0.02589 |
| mid | 35 | 29 | 6 | 0 | 0.02860 | 0.02855 |
| grad | 56 | 48 | 8 | 0 | 0.02134 | 0.02131 |
| big | 21 | 17 | 5 | 1 | 0.03045 | 0.03044 |
| mid, J5L | 39 | 31 | 8 | 0 | 0.02857 | 0.02853 |
<!-- /table:repair_accuracy -->

<!-- table:repair_cost -->
| sample | flagged | fits | changed | evals | per_fit_med | per_fit_mean | share | cap | verdict |
|---|---|---|---|---|---|---|---|---|---|
| held3 | 194 | 1214 | 45 | 842,911 | 662 | 694 | 3.59 | 2.5 | cost over the cap |
| held | 183 | 1026 | 41 | 719,418 | 662 | 701 | 3.14 | 2.5 | cost over the cap |
| mid | 142 | 635 | 37 | 451,286 | 662 | 711 | 4.26 | 2.5 | cost over the cap |
| grad | 238 | 965 | 32 | 640,885 | 662 | 664 | 7.09 | 3.5 | cost over the cap |
| big | 56 | 162 | 10 | 123,781 | 662 | 764 | 3.47 | 2.5 | cost over the cap |
| mid, J5L | 145 | 650 | 24 | 469,466 | 662 | 722 | 4.70 | n/a | measurement |
<!-- /table:repair_cost -->

<!-- table:repair_vs_s8 -->
| sample | j5_fixed | j5_broken | j5_net | r_fixed | r_broken | r_net |
|---|---|---|---|---|---|---|
| held3 | 24 | 17 | -7 | 25 | 15 | -10 |
| held | 26 | 23 | -3 | 29 | 17 | -12 |
| mid | 9 | 13 | 4 | 12 | 10 | -2 |
| grad | 31 | 24 | -7 | 36 | 21 | -15 |
| big | 8 | 5 | -3 | 13 | 6 | -7 |
<!-- /table:repair_vs_s8 -->
