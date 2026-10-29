<!-- table:s2_e0_evals -->
| arm | local_mean | global_mean | runtime_median | workers |
|---|---|---|---|---|
| key off | 4122.9 | 44499.9 | 2.44 | 6 |
| FusedCost 1 | 4119.2 | 44499.8 | 1.01 | 6 |
| BatchedQuickMC 1 | 4133.6 | 44501.1 | 1.60 | 6 |
| BatchedQuickMC 1, FusedCost 1 | 4132.7 | 44501.1 | 0.98 | 6 |
<!-- /table:s2_e0_evals -->

<!-- table:s2_e0_paired -->
| pair | n | wrong_a | wrong_b | fixed | broken | both_wrong | both_right | same_final |
|---|---|---|---|---|---|---|---|---|
| key off vs FusedCost 1 | 200 | 35 | 35 | 0 | 0 | 35 | 165 | 182 |
| BatchedQuickMC 1 vs BatchedQuickMC 1 + FusedCost 1 | 200 | 35 | 35 | 0 | 0 | 35 | 165 | 182 |
| key off vs BatchedQuickMC 1 (another draw, for scale) | 200 | 35 | 35 | 18 | 18 | 17 | 147 | 0 |
<!-- /table:s2_e0_paired -->

<!-- table:s2_e0_trace -->
| pair | voxel | change | reproduced | shadowed | differing | first_call | first_path | radius | fields | causes | within_rule | near_px |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| off vs on | 7336 | both wrong | True | 3165 | 87 | 55 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| off vs on | 16269 | both right | True | 4779 | 114 | 101 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0004 |
| off vs on | 20406 | both right | True | 4251 | 96 | 1 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 10531 | both right | True | 4958 | 126 | 30 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 16905 | both right | True | 4415 | 113 | 11 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0013 |
| off vs on | 14240 | both right | True | 5040 | 138 | 54 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| off vs on | 17977 | both right | True | 4070 | 102 | 117 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 10507 | both right | True | 3690 | 111 | 17 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0005 |
| off vs on | 22302 | both wrong | True | 3597 | 77 | 37 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 14497 | both wrong | True | 3461 | 109 | 11 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 12015 | both right | True | 4151 | 100 | 34 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 9525 | both right | True | 5674 | 145 | 30 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| off vs on | 15019 | both right | True | 5067 | 137 | 1 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| off vs on | 177 | both right | True | 4031 | 104 | 108 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| off vs on | 10864 | both right | True | 4519 | 108 | 30 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0006 |
| off vs on | 15620 | both right | True | 3439 | 77 | 17 | evaluate | 0 | pixel_on_detector, quality | pixel_integer | True | 0.0001 |
| off vs on | 22863 | both right | True | 4005 | 98 | 51 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| off vs on | 6596 | both right | True | 4527 | 144 | 17 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 2177 | both right | True | 3491 | 87 | 13 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 14217 | both wrong | True | 4247 | 112 | 35 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 5384 | both right | True | 3355 | 71 | 2 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 21741 | both right | True | 3653 | 92 | 138 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 17771 | both right | True | 3846 | 95 | 6 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 24339 | both right | True | 4944 | 113 | 110 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0009 |
| bqmc vs bqmc_on | 2807 | both right | True | 3747 | 106 | 21 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 1935 | both right | True | 4183 | 118 | 57 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 2735 | both right | True | 4714 | 99 | 86 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 22517 | both right | True | 3406 | 71 | 61 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 21726 | both right | True | 4147 | 119 | 90 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0015 |
| bqmc vs bqmc_on | 628 | both right | True | 3728 | 87 | 101 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 10007 | both wrong | True | 3587 | 95 | 38 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 24420 | both right | True | 3869 | 82 | 28 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0002 |
| bqmc vs bqmc_on | 20976 | both right | True | 3934 | 106 | 3 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 77 | both right | True | 4309 | 114 | 123 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
| bqmc vs bqmc_on | 18723 | both right | True | 4197 | 116 | 32 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0000 |
| bqmc vs bqmc_on | 14596 | both right | True | 4315 | 124 | 10 | evaluate | 0 | pixel_on_detector | pixel_integer | True | 0.0001 |
<!-- /table:s2_e0_trace -->

<!-- table:s2_e0_trace_pairs -->
| pair | voxels | changes | reproduced | shadowed | differing | within_rule | quality_differs | causes | max_near_px |
|---|---|---|---|---|---|---|---|---|---|
| bqmc vs bqmc_on | 18 | both right 16, both wrong 2 | 18 | 71672 | 1807 | 1807 | 198 | pixel_integer | 0.0018 |
| off vs on | 18 | both right 15, both wrong 3 | 18 | 76839 | 1986 | 1986 | 223 | pixel_integer | 0.0016 |
<!-- /table:s2_e0_trace_pairs -->

<!-- table:s2_e0_wrong -->
| arm | wrong | n | rate | wilson_low | wilson_high |
|---|---|---|---|---|---|
| key off | 35 | 200 | 0.175 | 0.129 | 0.234 |
| FusedCost 1 | 35 | 200 | 0.175 | 0.129 | 0.234 |
| BatchedQuickMC 1 | 35 | 200 | 0.175 | 0.129 | 0.234 |
| BatchedQuickMC 1, FusedCost 1 | 35 | 200 | 0.175 | 0.129 | 0.234 |
<!-- /table:s2_e0_wrong -->

<!-- table:s2_lib_batched -->
| set | K | B | compact | fused | fixed_gpu | gpu_over_fused | evaluate_torch | evaluate_fused |
|---|---|---|---|---|---|---|---|---|
| three0 | 112 | 1 | 251.2 | 16.9 | 2154.8 | 127.50 | 250.0 | 16.1 |
| three0 | 112 | 8 | 85.8 | 14.8 | 319.1 | 21.60 | 250.0 | 16.1 |
| three0 | 112 | 16 | 73.4 | 14.4 | 181.3 | 12.55 | 250.0 | 16.1 |
| three0 | 112 | 32 | 66.8 | 14.5 | 116.2 | 8.03 | 250.0 | 16.1 |
| three0 | 112 | 64 | 63.5 | 14.4 | 89.5 | 6.20 | 250.0 | 16.1 |
| three0 | 112 | 128 | 61.7 | 14.5 | 64.8 | 4.48 | 250.0 | 16.1 |
| three0 | 112 | 256 | 61.0 | 14.5 | 58.9 | 4.06 | 250.0 | 16.1 |
| three0 | 112 | 512 | 61.7 | 14.4 | 49.6 | 3.45 | 250.0 | 16.1 |
| three0_search | 868 | 1 | 618.3 | 95.9 | 2307.0 | 24.06 | 544.8 | 95.4 |
| three0_search | 868 | 8 | 448.8 | 93.1 | 435.3 | 4.68 | 544.8 | 95.4 |
| three0_search | 868 | 16 | 436.7 | 96.7 | 317.9 | 3.29 | 544.8 | 95.4 |
| three0_search | 868 | 32 | 433.1 | 98.1 | 271.6 | 2.77 | 544.8 | 95.4 |
| three0_search | 868 | 64 | 431.3 | 98.1 | 236.6 | 2.41 | 544.8 | 95.4 |
| three0_search | 868 | 128 | 431.9 | 96.2 | 227.2 | 2.36 | 544.8 | 95.4 |
| three0_search | 868 | 256 | 426.9 | 94.3 | 227.8 | 2.42 | 544.8 | 95.4 |
| three0_search | 868 | 512 | 421.0 | 90.5 | 229.8 | 2.54 | 544.8 | 95.4 |
| e0_0 | 112 | 1 | 246.5 | 17.6 | 2228.9 | 126.36 | 248.3 | 16.9 |
| e0_0 | 112 | 8 | 82.2 | 15.5 | 324.4 | 20.94 | 248.3 | 16.9 |
| e0_0 | 112 | 16 | 69.0 | 15.3 | 185.7 | 12.12 | 248.3 | 16.9 |
| e0_0 | 112 | 32 | 62.6 | 15.2 | 120.4 | 7.94 | 248.3 | 16.9 |
| e0_0 | 112 | 64 | 59.0 | 15.2 | 80.2 | 5.28 | 248.3 | 16.9 |
| e0_0 | 112 | 128 | 57.1 | 15.1 | 64.7 | 4.28 | 248.3 | 16.9 |
| e0_0 | 112 | 256 | 56.4 | 15.0 | 59.6 | 3.97 | 248.3 | 16.9 |
| e0_0 | 112 | 512 | 56.2 | 15.1 | 48.2 | 3.20 | 248.3 | 16.9 |
<!-- /table:s2_lib_batched -->

<!-- table:s2_lib_intermediates -->
| dataset | intermediate | n | equal | equal_frac | max_ulp | max_abs |
|---|---|---|---|---|---|---|
| three0 | g_lab | 336000 | 312648 | 0.930 | 1543 | 4.77e-07 |
| three0 | omega_bragg | 221820 | 206137 | 0.929 | 1024 | 4.77e-07 |
| three0 | omega | 202504 | 187762 | 0.927 | 1024 | 4.77e-07 |
| three0 | eta | 202504 | 169248 | 0.836 | 33 | 2.38e-07 |
| three0 | pixel | 1221696 | 1035205 | 0.847 | 254 | 3.88e-03 |
| three0 | t | 610848 | 580405 | 0.950 | 5 | 1.91e-06 |
<!-- /table:s2_lib_intermediates -->

<!-- table:s2_lib_parity -->
| dataset | modes | n | differing | silent | flipped_coords | failing | causes | max_near_px | max_margin_ulp | several_coords | frame_differ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| three0 | triangle (5 modes) | 1000 | 16 | 25-27 | 41-43 | 0 | pixel_integer | 0.0009 | 148 | 0 | 0 |
| three0 | screen (4 modes) | 1000 | 0 | 1-2 | 1-2 | 0 | pixel_integer | 0.0000 | 1 | 0 | 0 |
| full15901 | triangle (5 modes) | 1000 | 28 | 29 | 60 | 0 | pixel_integer | 0.0009 | 49 | 0 | 0 |
| full15901 | screen (4 modes) | 1000 | 0 | 4-6 | 4-6 | 0 | pixel_integer | 0.0004 | 20 | 0 | 0 |
| full17034 | triangle (5 modes) | 1000 | 29 | 33 | 64 | 0 | pixel_integer | 0.0014 | 265 | 0 | 0 |
| full17034 | screen (4 modes) | 1000 | 0 | 2-5 | 2-5 | 0 | pixel_integer | 0.0010 | 76 | 0 | 0 |
| full18558 | triangle (5 modes) | 1000 | 19-20 | 32 | 54-55 | 0 | pixel_integer | 0.0006 | 312 | 0 | 0 |
| full18558 | screen (4 modes) | 1000 | 0 | 4-5 | 4-5 | 0 | pixel_integer | 0.0007 | 88 | 0 | 0 |
<!-- /table:s2_lib_parity -->

<!-- table:s2_lib_parity_local -->
| path | n | differing | wilson_low | wilson_high | pixel_on_detector | pixel_overlap | quality | peak_overlap | peak_on_detector | flipped_coords | nearer_torch | nearer_fused | side_torch | side_fused |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| library (plain g_lab) | 4000 | 93 | 0.019 | 0.028 | 92 | 5 | 5 | 0 | 0 | 222 | 113 | 109 | 109 | 113 |
<!-- /table:s2_lib_parity_local -->

<!-- table:s2_lib_per_call -->
| arm | torch | fused | fused_range | c_only | python_share | speedup | step0_prototype |
|---|---|---|---|---|---|---|---|
| triangle_q8_default | 243.7 | 17.2 | 16.8-18.5 | 13.7 | 3.5 | 14.2 | 14.4 |
| triangle_q8_cpp | 242.6 | 17.0 | 16.7-17.1 | 13.6 | 3.4 | 14.3 | 14.2 |
| screen_r3_q5_default | 206.6 | 4.6 | 4.4-4.7 | 2.4 | 2.1 | 45.2 | 3.3 |
| screen_r3_q5_cpp | 206.2 | 4.4 | 4.4-4.5 | 2.4 | 2.0 | 46.9 | 3.3 |
<!-- /table:s2_lib_per_call -->

<!-- table:s2_lib_search_e0 -->
| arm | wall | vs_default | evaluate | evaluate_n | us_per_eval | quick_mc | quick_mc_n | quick_mc_us | screen | wrong |
|---|---|---|---|---|---|---|---|---|---|---|
| default keys | 16.49 | 1.00 | 10.18 | 39759 | 256 | 7.04 | 27636 | 255 | 5.31 | 2 of 10 |
| BatchedQuickMC 1 | 10.91 | 1.51 | 2.89 | 11219 | 258 | 1.79 | 27576 | 65 | 5.40 | 2 of 10 |
| FusedCost 1 | 7.02 | 2.35 | 0.84 | 39869 | 21 | 0.58 | 27636 | 21 | 5.30 | 2 of 10 |
| BatchedQuickMC 1, FusedCost 1 | 6.91 | 2.39 | 0.23 | 11219 | 21 | 0.58 | 27576 | 21 | 5.32 | 2 of 10 |
<!-- /table:s2_lib_search_e0 -->

<!-- table:s2_lib_search_three -->
| arm | wall | vs_default | evaluate | evaluate_n | us_per_eval | quick_mc | quick_mc_n | quick_mc_us | screen | wrong |
|---|---|---|---|---|---|---|---|---|---|---|
| default keys | 5.50 | 1.00 | 4.67 | 8311 | 562 | 3.63 | 6456 | 562 | 0.64 | 1 of 1 |
| BatchedQuickMC 1 | 4.60 | 1.20 | 1.02 | 1812 | 563 | 2.79 | 6408 | 435 | 0.64 | 0 of 1 |
| FusedCost 1 | 1.77 | 3.11 | 0.96 | 8311 | 116 | 0.75 | 6456 | 117 | 0.65 | 1 of 1 |
| BatchedQuickMC 1, FusedCost 1 | 1.81 | 3.04 | 0.21 | 1812 | 114 | 0.82 | 6408 | 129 | 0.65 | 0 of 1 |
<!-- /table:s2_lib_search_three -->
