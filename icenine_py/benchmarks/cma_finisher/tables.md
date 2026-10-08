<!-- table:e0_accuracy -->
| variant | method | n | wrong | med_err | med_both | cost |
|---|---|---|---|---|---|---|
| clean | mc | 200 | 67/200 = 0.335 [0.273, 0.403] | 0.0301 | 0.0301 | 0.3542 |
| clean | cma | 200 | 66/200 = 0.330 [0.269, 0.398] | 0.0019 | 0.0019 | 0.2959 |
| realistic | mc | 200 | 52/200 = 0.260 [0.204, 0.325] | 0.0278 | 0.0278 | 0.4259 |
| realistic | cma | 200 | 44/200 = 0.220 [0.168, 0.282] | 0.0021 | 0.0021 | 0.3693 |
<!-- /table:e0_accuracy -->

<!-- table:e0_evals_wall -->
| variant | method | glob | loc | find | wall | wall_med |
|---|---|---|---|---|---|---|
| clean | mc | 44478 | 4711 | 2.54 | 24.72 | 24.54 |
| clean | cma | 44478 | 4805 | 1.90 | 24.90 | 24.85 |
| realistic | mc | 44778 | 7628 | 4.82 | 26.81 | 26.59 |
| realistic | cma | 44778 | 9461 | 4.82 | 28.24 | 27.60 |
<!-- /table:e0_evals_wall -->

<!-- table:timing_single_worker -->
| variant | n | mc_mean_s | cma_mean_s | ratio_of_means | median_paired_ratio | mc_mean_evals | cma_mean_evals |
|---|---|---|---|---|---|---|---|
| clean | 20 | 19.47 | 19.81 | 1.017 | 1.033 | 49083 | 49619 |
| realistic | 20 | 21.57 | 22.60 | 1.048 | 1.040 | 52801 | 54669 |
<!-- /table:timing_single_worker -->

<!-- table:validate_lib -->
| method | n | median error (deg) | fraction < 0.02 deg (Wilson 95%) | median evals (local_optimization rows) |
|---|---|---|---|---|
| start (net x3) | 50 | 0.0690 | 12/50 [0.14, 0.37] | n/a |
| B3 cma_02, 250 evals | 50 | 0.0034 | 45/50 [0.79, 0.96] | n/a |
| library CMAOptimizer, 250 evals | 50 | 0.0034 | 45/50 [0.79, 0.96] | n/a |
| B3 cma_02, 1000 evals | 50 | 0.0025 | 48/50 [0.87, 0.99] | n/a |
| library CMAOptimizer, 1000 evals | 50 | 0.0025 | 48/50 [0.87, 0.99] | n/a |
| refine_from_candidates, local_optimizer cma (1000 + final) | 50 | 0.0035 | 45/50 [0.79, 0.96] | n/a |
| local_optimization, mc (VarianceMinimizing) | 50 | 0.0690 | 12/50 [0.14, 0.37] | 662 |
| local_optimization, cma | 50 | 0.0028 | 48/50 [0.87, 0.99] | 1001 |
<!-- /table:validate_lib -->
