<!-- table:e0_accuracy -->
| variant | method | n | wrong | med_err | med_both | cost |
|---|---|---|---|---|---|---|
| clean | mc | 200 | 40/200 = 0.200 [0.150, 0.261] | 0.0107 | 0.0107 | 0.2340 |
| clean | cma | 200 | 40/200 = 0.200 [0.150, 0.261] | 0.0018 | 0.0018 | 0.2003 |
| realistic | mc | 200 | 25/200 = 0.125 [0.086, 0.178] | 0.0087 | 0.0087 | 0.3358 |
| realistic | cma | 200 | 24/200 = 0.120 [0.082, 0.172] | 0.0021 | 0.0021 | 0.3196 |
<!-- /table:e0_accuracy -->

<!-- table:e0_evals_wall -->
| variant | method | glob | loc | find | wall | wall_med |
|---|---|---|---|---|---|---|
| clean | mc | 44479 | 3963 | 2.40 | 67.72 | 76.94 |
| clean | cma | 44479 | 4892 | 1.74 | 69.74 | 78.28 |
| realistic | mc | 44780 | 6278 | 4.68 | 72.73 | 82.40 |
| realistic | cma | 44780 | 9711 | 4.68 | 79.16 | 88.22 |
<!-- /table:e0_evals_wall -->

<!-- table:timing_single_worker -->
| variant | n | mc_mean_s | cma_mean_s | ratio_of_means | median_paired_ratio | mc_mean_evals | cma_mean_evals |
|---|---|---|---|---|---|---|---|
| clean | 20 | 19.47 | 19.81 | 1.017 | 1.033 | 49083 | 49619 |
| realistic | 20 | 21.57 | 22.60 | 1.048 | 1.040 | 52801 | 54669 |
<!-- /table:timing_single_worker -->
