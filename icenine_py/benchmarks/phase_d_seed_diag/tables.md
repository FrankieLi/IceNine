<!-- table:seed_diag_cpp -->
| config | images | voxel | n_runs | adap_s_mean | adap_evals_mean | us_per_eval | variance_steps_median | variance_steps_max |
|---|---|---|---|---|---|---|---|---|
| Phase D search keys (MaxLocalResolution 3, 200 steps, 2 restarts, 30 candidates) | full | 15901 | 16 | 13.3 | 258,996 | 51.3 | 1.3e+03 | 2130 |
| Phase D search keys (MaxLocalResolution 3, 200 steps, 2 restarts, 30 candidates) | full | 22867 | 16 | 12.8 | 262,537 | 48.8 | 1.28e+03 | 1930 |
| C++ example config (MaxLocalResolution 5, 300 steps, 3 restarts, 50 candidates) | full | 15901 | 16 | 13.9 | 282,813 | 49.0 | 2.19e+04 | 28840 |
| C++ example config (MaxLocalResolution 5, 300 steps, 3 restarts, 50 candidates) | full | 22867 | 16 | 14.1 | 286,758 | 49.1 | 1.95e+04 | 33980 |
| Phase D search keys, isolated one-voxel images | isolated | 15901 | 16 | 0.6 | 47,833 | 12.5 | 325 | 460 |
<!-- /table:seed_diag_cpp -->

<!-- table:seed_diag_isolated_vs_full -->
| voxel | iso_evals | iso_wall_s | iso_us | iso_err | iso_var_evals | full_evals | full_wall_s | full_us | full_err | evals_ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 50,091 | 22 | 434 | 0.129 | 1,387 | 1,687,675 | 897 | 503 | 0.042 | 33.7 |
| 17034 | 48,397 | 20 | 401 | 59.981 | 452 | 952,819 | 509 | 505 | 0.030 | 19.7 |
| 18558 | 48,096 | 19 | 396 | 59.977 | 485 | 463,626 | 245 | 499 | 0.006 | 9.6 |
| 19889 | 49,292 | 20 | 401 | 0.003 | 1,453 | 626,184 | 333 | 502 | 0.024 | 12.7 |
| 22867 | 48,340 | 20 | 405 | 25.267 | 298 | 1,339,563 | 753 | 533 | 0.024 | 27.7 |
| 5242 | 49,174 | 20 | 404 | 0.049 | 870 | 1,063,389 | 589 | 524 | 0.028 | 21.6 |
<!-- /table:seed_diag_isolated_vs_full -->

<!-- table:seed_diag_levels -->
| level | n_fz | n_returned | n_kept | discrete_evals | n_recip | diameter_deg |
|---|---|---|---|---|---|---|
| 0 | 4,886 | 10,267 | 2,566 | 78,675 | 26 | 5.00 |
| 1 | 2,566 | 3,452 | 862 | 45,191 | 50 | 3.33 |
| 2 | 862 | 1,080 | 270 | 15,397 | 64 | 2.22 |
| 3 | 270 | 336 | 84 | 4,853 | 112 | 1.48 |
<!-- /table:seed_diag_levels -->

<!-- table:seed_diag_options_interior -->
| option | n | wall_s | evals | err_med | err_max | wrong |
|---|---|---|---|---|---|---|
| base (mc, as Phase D config) | 6 | 554 | 1,022,209 | 0.026 | 0.042 | 0/6 [0.00, 0.39] |
| variance cap 2000 steps | 6 | 169 | 316,678 | 0.029 | 0.061 | 0/6 [0.00, 0.39] |
| CMA-ES finisher (LocalOptimizer cma) | 6 | 178 | 326,268 | 0.025 | 0.046 | 0/6 [0.00, 0.39] |
| cap 2000 + top 3000 at level 0 + keep 1/8 | 6 | 68 | 129,218 | 0.029 | 0.053 | 0/6 [0.00, 0.39] |
<!-- /table:seed_diag_options_interior -->

<!-- table:seed_diag_options_random -->
| option | n | wall_s | wall_max | evals | err_med | err_max | wrong |
|---|---|---|---|---|---|---|---|
| variance cap 2000 | 18 | 172 | 182 | 317,146 | 0.028 | 52.342 | 1/18 [0.01, 0.26] |
| lean (cap 2000, top 3000, keep 1/8) | 18 | 67 | 71 | 128,857 | 0.028 | 59.940 | 2/18 [0.03, 0.33] |
<!-- /table:seed_diag_options_random -->

<!-- table:seed_diag_ranks -->
| voxel | n | n_near | best_near_rank_discrete | best_near_rank_after_quickmc | n_near_in_kept_quarter |
|---|---|---|---|---|---|
| 15901 | 10,190 | 4 | 0 | 0 | 4 |
| 17034 | 10,382 | 3 | 27 | 68 | 3 |
| 18558 | 10,271 | 4 | 9 | 0 | 3 |
| 19889 | 10,308 | 2 | 7 | 0 | 2 |
| 22867 | 10,236 | 3 | 1213 | 74 | 3 |
| 5242 | 10,213 | 3 | 507 | 35 | 3 |
<!-- /table:seed_diag_ranks -->

<!-- table:seed_diag_run_estimates -->
| option | per_seed_s | seeds_low | seeds_mid | seeds_high | hours_low | hours_mid | hours_high |
|---|---|---|---|---|---|---|---|
| base (as Phase D config) | 486 | 350 | 395 | 504 | 47.2 | 53.3 | 68.0 |
| variance cap 2000 | 150 | 350 | 395 | 504 | 14.6 | 16.5 | 21.1 |
| CMA-ES finisher (no cap needed) | 157 | 350 | 395 | 504 | 15.2 | 17.2 | 21.9 |
| lean: cap 2000 + top 3000 at level 0 + keep 1/8 | 62 | 350 | 395 | 504 | 6.0 | 6.8 | 8.6 |
<!-- /table:seed_diag_run_estimates -->

<!-- table:seed_diag_screen -->
| images | voxel | screened | passed | pass_rate | n_returned |
|---|---|---|---|---|---|
| full | 15901 | 43,974 | 35,267 | 0.802 | 10,190 |
| full | 17034 | 43,974 | 35,780 | 0.814 | 10,382 |
| full | 18558 | 43,974 | 33,617 | 0.764 | 10,271 |
| full | 19889 | 43,974 | 34,689 | 0.789 | 10,308 |
| full | 22867 | 43,974 | 35,512 | 0.808 | 10,236 |
| full | 5242 | 43,974 | 33,342 | 0.758 | 10,213 |
| isolated | 15901 | 43,974 | 224 | 0.005 | 207 |
| isolated | 17034 | 43,974 | 193 | 0.004 | 182 |
| isolated | 18558 | 43,974 | 169 | 0.004 | 166 |
| isolated | 19889 | 43,974 | 164 | 0.004 | 161 |
| isolated | 22867 | 43,974 | 183 | 0.004 | 172 |
| isolated | 5242 | 43,974 | 201 | 0.005 | 189 |
<!-- /table:seed_diag_screen -->

<!-- table:seed_diag_seed_count -->
| p_fail | mean | lo | hi |
|---|---|---|---|
| 0.0 | 356 | 348 | 367 |
| 0.05 | 377 | 366 | 394 |
| 0.1 | 395 | 381 | 417 |
| 0.2 | 443 | 426 | 464 |
| 0.3 | 504 | 474 | 532 |
<!-- /table:seed_diag_seed_count -->

<!-- table:seed_diag_stage_breakdown -->
| voxel | wall_s | evals | disc_L0 | quick_L0 | rest_L1_3 | find | variance | var_share | us | find_n | conv | err | q_true |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 897 | 1,687,675 | 79,241 | 112,090 | 118,226 | 201 | 1,377,916 | 0.82 | 503 | 1/344 | True | 0.042 | 0.947 |
| 17034 | 509 | 952,819 | 79,754 | 114,202 | 120,721 | 5,849 | 632,292 | 0.66 | 505 | 30/341 | False | 0.030 | 0.922 |
| 18558 | 245 | 463,626 | 77,591 | 112,981 | 118,251 | 603 | 154,199 | 0.33 | 499 | 3/330 | True | 0.006 | 0.937 |
| 19889 | 333 | 626,184 | 78,663 | 113,388 | 119,796 | 5,466 | 308,870 | 0.49 | 502 | 30/334 | False | 0.024 | 0.920 |
| 22867 | 753 | 1,339,563 | 79,486 | 112,596 | 119,087 | 5,579 | 1,022,814 | 0.76 | 533 | 30/337 | False | 0.024 | 0.941 |
| 5242 | 589 | 1,063,389 | 77,316 | 112,343 | 117,868 | 5,561 | 750,300 | 0.71 | 524 | 30/333 | False | 0.028 | 0.908 |
| mean | 554 | 1,022,209 | 78,675 | 112,933 | 118,992 | 3,876 | 707,732 | 0.63 | 511 |  |  | n/a | n/a |
<!-- /table:seed_diag_stage_breakdown -->

<!-- table:seed_diag_variance_trace -->
| voxel | images | box_deg | ended | n_runs | p_low | var_median |
|---|---|---|---|---|---|---|
| 15901 | full | 0.33 | cap | 3,000 | 0.0003 | 0.0213 |
| 15901 | full | 0.22 | cap | 3,000 | 0.0003 | 0.0175 |
| 15901 | full | 0.20 | cap | 3,000 | 0.0003 | 0.0156 |
| 15901 | isolated | 0.33 | terminated | 189 | 0.1058 | 0.0069 |
| 15901 | isolated | 0.22 | terminated | 1,359 | 0.0147 | 0.0134 |
| 15901 | isolated | 0.20 | terminated | 2,527 | 0.0079 | 0.0149 |
<!-- /table:seed_diag_variance_trace -->
