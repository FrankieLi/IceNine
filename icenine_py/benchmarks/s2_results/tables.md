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
| three0 | 112 | 1 | 245.7 | 16.6 | 2156.8 | 129.68 | 244.1 | 15.9 |
| three0 | 112 | 8 | 87.8 | 14.6 | 317.3 | 21.75 | 244.1 | 15.9 |
| three0 | 112 | 16 | 74.6 | 14.4 | 184.7 | 12.83 | 244.1 | 15.9 |
| three0 | 112 | 32 | 68.4 | 14.3 | 119.4 | 8.36 | 244.1 | 15.9 |
| three0 | 112 | 64 | 65.0 | 14.2 | 97.3 | 6.83 | 244.1 | 15.9 |
| three0 | 112 | 128 | 63.4 | 14.2 | 65.5 | 4.60 | 244.1 | 15.9 |
| three0 | 112 | 256 | 62.7 | 14.2 | 59.2 | 4.17 | 244.1 | 15.9 |
| three0 | 112 | 512 | 63.4 | 14.1 | 48.4 | 3.44 | 244.1 | 15.9 |
| three0_search | 868 | 1 | 628.2 | 93.3 | 2303.5 | 24.69 | 543.8 | 95.5 |
| three0_search | 868 | 8 | 459.2 | 92.7 | 425.0 | 4.59 | 543.8 | 95.5 |
| three0_search | 868 | 16 | 446.6 | 95.3 | 316.5 | 3.32 | 543.8 | 95.5 |
| three0_search | 868 | 32 | 442.6 | 97.2 | 272.7 | 2.80 | 543.8 | 95.5 |
| three0_search | 868 | 64 | 441.3 | 97.3 | 235.6 | 2.42 | 543.8 | 95.5 |
| three0_search | 868 | 128 | 440.5 | 95.5 | 228.0 | 2.39 | 543.8 | 95.5 |
| three0_search | 868 | 256 | 437.5 | 95.3 | 228.0 | 2.39 | 543.8 | 95.5 |
| three0_search | 868 | 512 | 431.5 | 92.8 | 228.4 | 2.46 | 543.8 | 95.5 |
| e0_0 | 112 | 1 | 248.6 | 17.6 | 2208.6 | 125.20 | 248.0 | 16.9 |
| e0_0 | 112 | 8 | 83.3 | 15.4 | 318.5 | 20.64 | 248.0 | 16.9 |
| e0_0 | 112 | 16 | 70.4 | 15.3 | 187.8 | 12.30 | 248.0 | 16.9 |
| e0_0 | 112 | 32 | 63.6 | 15.2 | 117.9 | 7.75 | 248.0 | 16.9 |
| e0_0 | 112 | 64 | 60.4 | 15.1 | 78.8 | 5.22 | 248.0 | 16.9 |
| e0_0 | 112 | 128 | 58.9 | 15.0 | 70.4 | 4.69 | 248.0 | 16.9 |
| e0_0 | 112 | 256 | 57.7 | 15.1 | 56.9 | 3.77 | 248.0 | 16.9 |
| e0_0 | 112 | 512 | 57.9 | 15.1 | 53.2 | 3.51 | 248.0 | 16.9 |
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
| triangle_q8_default | 235.7 | 15.9 | 15.4-16.1 | 13.8 | 2.0 | 14.9 | 14.4 |
| triangle_q8_cpp | 235.0 | 15.6 | 15.3-15.8 | 13.7 | 1.9 | 15.1 | 14.2 |
| screen_r3_q5_default | 202.1 | 4.4 | 4.3-4.5 | 2.5 | 1.9 | 45.5 | 3.3 |
| screen_r3_q5_cpp | 202.1 | 4.3 | 4.2-4.4 | 2.4 | 1.9 | 46.8 | 3.3 |
<!-- /table:s2_lib_per_call -->

<!-- table:s2_lib_search_e0 -->
| arm | wall | vs_default | evaluate | evaluate_n | us_per_eval | quick_mc | quick_mc_n | quick_mc_us | screen | wrong |
|---|---|---|---|---|---|---|---|---|---|---|
| default keys | 16.10 | 1.00 | 9.95 | 39759 | 250 | 6.91 | 27636 | 250 | 5.16 | 2 of 10 |
| BatchedQuickMC 1 | 10.52 | 1.53 | 2.82 | 11219 | 251 | 1.73 | 27576 | 63 | 5.16 | 2 of 10 |
| FusedCost 1 | 6.84 | 2.35 | 0.80 | 39869 | 20 | 0.56 | 27636 | 20 | 5.16 | 2 of 10 |
| BatchedQuickMC 1, FusedCost 1 | 6.72 | 2.39 | 0.22 | 11219 | 20 | 0.56 | 27576 | 20 | 5.16 | 2 of 10 |
<!-- /table:s2_lib_search_e0 -->

<!-- table:s2_lib_search_three -->
| arm | wall | vs_default | evaluate | evaluate_n | us_per_eval | quick_mc | quick_mc_n | quick_mc_us | screen | wrong |
|---|---|---|---|---|---|---|---|---|---|---|
| default keys | 5.40 | 1.00 | 4.60 | 8311 | 553 | 3.58 | 6456 | 554 | 0.62 | 1 of 1 |
| BatchedQuickMC 1 | 4.50 | 1.20 | 1.00 | 1812 | 550 | 2.72 | 6408 | 425 | 0.62 | 0 of 1 |
| FusedCost 1 | 1.72 | 3.14 | 0.94 | 8311 | 113 | 0.74 | 6456 | 114 | 0.63 | 1 of 1 |
| BatchedQuickMC 1, FusedCost 1 | 1.76 | 3.07 | 0.20 | 1812 | 110 | 0.80 | 6408 | 125 | 0.63 | 0 of 1 |
<!-- /table:s2_lib_search_three -->
