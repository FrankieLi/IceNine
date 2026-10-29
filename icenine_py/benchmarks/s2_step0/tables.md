<!-- table:s2_step0_callers_e0 -->
| arm | caller | calls | seconds | share | us_per_call |
|---|---|---|---|---|---|
| default keys | screen (BatchedScreen, not in v1) | 0 | 5.16 | 0.318 | n/a |
| default keys | recost | 3711 | 0.96 | 0.059 | 259 |
| default keys | quick_mc | 27636 | 7.00 | 0.432 | 253 |
| default keys | find_optimal | 4938 | 1.23 | 0.076 | 250 |
| default keys | variance | 3464 | 0.87 | 0.053 | 250 |
| default keys | final | 10 | 0.00 | 0.000 | 253 |
| default keys | batched quick-MC cost | 0 | 0.00 | 0.000 | n/a |
| default keys | covered evaluate (all but screen) | 39759 | 10.06 | 0.621 | 253 |
| BatchedQuickMC 1 | screen (BatchedScreen, not in v1) | 0 | 5.15 | 0.490 | n/a |
| BatchedQuickMC 1 | recost | 3650 | 0.95 | 0.090 | 260 |
| BatchedQuickMC 1 | quick_mc | 0 | 0.00 | 0.000 | n/a |
| BatchedQuickMC 1 | find_optimal | 3996 | 1.00 | 0.095 | 251 |
| BatchedQuickMC 1 | variance | 3563 | 0.89 | 0.085 | 250 |
| BatchedQuickMC 1 | final | 10 | 0.00 | 0.000 | 254 |
| BatchedQuickMC 1 | batched quick-MC cost | 27576 | 1.71 | 0.163 | 62 |
| BatchedQuickMC 1 | covered evaluate (all but screen) | 11219 | 2.84 | 0.270 | 254 |
<!-- /table:s2_step0_callers_e0 -->

<!-- table:s2_step0_callers_three -->
| arm | caller | calls | seconds | share | us_per_call |
|---|---|---|---|---|---|
| default keys | screen (BatchedScreen, not in v1) | 0 | 0.63 | 0.115 | n/a |
| default keys | recost | 893 | 0.50 | 0.093 | 564 |
| default keys | quick_mc | 6456 | 3.62 | 0.665 | 560 |
| default keys | find_optimal | 311 | 0.16 | 0.030 | 528 |
| default keys | variance | 650 | 0.34 | 0.063 | 530 |
| default keys | final | 1 | 0.00 | 0.000 | 535 |
| default keys | batched quick-MC cost | 0 | 0.00 | 0.000 | n/a |
| default keys | covered evaluate (all but screen) | 8311 | 4.63 | 0.852 | 557 |
| BatchedQuickMC 1 | screen (BatchedScreen, not in v1) | 0 | 0.63 | 0.138 | n/a |
| BatchedQuickMC 1 | recost | 916 | 0.52 | 0.114 | 565 |
| BatchedQuickMC 1 | quick_mc | 0 | 0.00 | 0.000 | n/a |
| BatchedQuickMC 1 | find_optimal | 311 | 0.16 | 0.036 | 529 |
| BatchedQuickMC 1 | variance | 584 | 0.31 | 0.068 | 530 |
| BatchedQuickMC 1 | final | 1 | 0.00 | 0.000 | 528 |
| BatchedQuickMC 1 | batched quick-MC cost | 6408 | 2.76 | 0.610 | 431 |
| BatchedQuickMC 1 | covered evaluate (all but screen) | 1812 | 0.99 | 0.219 | 548 |
<!-- /table:s2_step0_callers_three -->

<!-- table:s2_step0_float64_error -->
| dataset | variant | torch_mean_px | fused_mean_px | torch_max_px | fused_max_px |
|---|---|---|---|---|---|
| three0 | portable | 2.87e-04 | 2.85e-04 | 0.0112 | 0.0112 |
| three0 | accelerate | 2.87e-04 | 2.87e-04 | 0.0112 | 0.0112 |
| full15901 | portable | 2.86e-04 | 2.84e-04 | 0.0058 | 0.0066 |
| full15901 | accelerate | 2.86e-04 | 2.86e-04 | 0.0058 | 0.0058 |
| full17034 | portable | 2.85e-04 | 2.84e-04 | 0.0066 | 0.0066 |
| full17034 | accelerate | 2.85e-04 | 2.86e-04 | 0.0066 | 0.0066 |
| full18558 | portable | 2.98e-04 | 2.95e-04 | 0.0669 | 0.0669 |
| full18558 | accelerate | 2.98e-04 | 2.99e-04 | 0.0669 | 0.0669 |
<!-- /table:s2_step0_float64_error -->

<!-- table:s2_step0_fmaf_probe -->
| product | order | n | equal | equal_frac | max_ulp |
|---|---|---|---|---|---|
| g_lab (R @ g_hkl.T).T | plain | 6720 | 6244 | 0.929 | 73 |
| g_lab (R @ g_hkl.T).T | fma_chain | 6720 | 5404 | 0.804 | 128 |
| g_lab (R @ g_hkl.T).T | fma_rev | 6720 | 4546 | 0.676 | 286 |
| torch.norm(g_lab) | plain | 2218 | 1926 | 0.868 | 1 |
| torch.norm(g_lab) | fma_chain | 2218 | 2218 | 1.000 | 0 |
| torch.norm(g_lab) | fma_rev | 2218 | 1802 | 0.812 | 1 |
| Rz @ base | plain | 39924 | 39924 | 1.000 | 0 |
| Rz @ base | fma_chain | 39924 | 39924 | 1.000 | 0 |
| Rz @ base | fma_rev | 39924 | 39924 | 1.000 | 0 |
| bmm(full_rot, normals) | plain | 13308 | 10007 | 0.752 | 8 |
| bmm(full_rot, normals) | fma_chain | 13308 | 13308 | 1.000 | 0 |
| bmm(full_rot, normals) | fma_rev | 13308 | 9280 | 0.697 | 9 |
| torch.norm(reflected) | plain | 4436 | 4374 | 0.986 | 1 |
| torch.norm(reflected) | fma_chain | 4436 | 4436 | 1.000 | 0 |
| torch.norm(reflected) | fma_rev | 4436 | 3721 | 0.839 | 1 |
| einsum vertices (4 x 4) | plain | 39924 | 39025 | 0.977 | 480 |
| einsum vertices (4 x 4) | fma_chain | 39924 | 39924 | 1.000 | 0 |
| einsum vertices (4 x 4) | fma_rev | 39924 | 37397 | 0.937 | 973 |
<!-- /table:s2_step0_fmaf_probe -->

<!-- table:s2_step0_intermediates -->
| variant | intermediate | n | equal | equal_frac | max_ulp | max_abs |
|---|---|---|---|---|---|---|
| portable | g_lab | 336000 | 270828 | 0.806 | 10201 | 9.54e-07 |
| portable | omega_bragg | 221820 | 197236 | 0.889 | 8704 | 7.01e-05 |
| portable | omega | 202504 | 179463 | 0.886 | 8704 | 7.01e-05 |
| portable | eta | 202504 | 152736 | 0.754 | 6648 | 6.51e-06 |
| portable | pixel | 1221696 | 931346 | 0.762 | 253 | 9.77e-03 |
| portable | t | 610848 | 562887 | 0.921 | 7 | 2.86e-06 |
| accelerate | g_lab | 336000 | 336000 | 1.000 | 0 | 0.00e+00 |
| accelerate | omega_bragg | 221820 | 212958 | 0.960 | 1024 | 4.77e-07 |
| accelerate | omega | 202504 | 194325 | 0.960 | 1024 | 4.77e-07 |
| accelerate | eta | 202504 | 178256 | 0.880 | 32 | 1.79e-07 |
| accelerate | pixel | 1221696 | 1121235 | 0.918 | 254 | 3.88e-03 |
| accelerate | t | 610848 | 593251 | 0.971 | 5 | 1.91e-06 |
<!-- /table:s2_step0_intermediates -->

<!-- table:s2_step0_parity -->
| dataset | variant | modes | differing | beyond_8ulp | no_peak | causes | max_margin_ulp | max_margin_px |
|---|---|---|---|---|---|---|---|---|
| three0 | portable | triangle (5 modes) | 25 | 5 | 0 | pixel_integer | 54 | 0.0016 |
| three0 | portable | screen (4 modes) | 0-1 | 0 | 0 | pixel_integer | 1 | 0.0001 |
| three0 | accelerate | triangle (5 modes) | 8 | 3 | 0 | pixel_integer | 43 | 0.0009 |
| three0 | accelerate | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
| full15901 | portable | triangle (5 modes) | 36-42 | 4-7 | 0 | pixel_integer | 49 | 0.0029 |
| full15901 | portable | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
| full15901 | accelerate | triangle (5 modes) | 15 | 3 | 0 | pixel_integer | 49 | 0.0015 |
| full15901 | accelerate | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
| full17034 | portable | triangle (5 modes) | 32-34 | 10-11 | 0 | pixel_integer | 197 | 0.0021 |
| full17034 | portable | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
| full17034 | accelerate | triangle (5 modes) | 12 | 5 | 0 | pixel_integer | 197 | 0.0020 |
| full17034 | accelerate | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
| full18558 | portable | triangle (5 modes) | 22-28 | 3-11 | 0 | pixel_integer | 312 | 0.0033 |
| full18558 | portable | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
| full18558 | accelerate | triangle (5 modes) | 9 | 2 | 0 | pixel_integer | 312 | 0.0012 |
| full18558 | accelerate | screen (4 modes) | 0 | 0 | 0 | - | n/a | n/a |
<!-- /table:s2_step0_parity -->

<!-- table:s2_step0_per_call -->
| arm | torch | fused | fused_range | fused_acc | c_only | speedup |
|---|---|---|---|---|---|---|
| triangle_q8_default | 234.1 | 14.4 | 13.9-14.6 | 14.4 | 13.4 | 16.2 |
| triangle_q8_cpp | 233.7 | 14.2 | 13.9-14.6 | 14.3 | 13.5 | 16.4 |
| screen_r3_q5_default | 202.0 | 3.3 | 3.1-3.4 | 3.3 | 2.5 | 61.4 |
| screen_r3_q5_cpp | 201.5 | 3.3 | 3.2-3.4 | 3.3 | 2.6 | 60.6 |
<!-- /table:s2_step0_per_call -->

<!-- table:s2_step0_search_e0 -->
| arm | wall | vs_default | evaluate | evaluate_n | us_per_eval | screen | bqmc | wrong |
|---|---|---|---|---|---|---|---|---|
| default keys | 16.20 | 1.00 | 10.06 | 39759 | 253 | 5.16 | 0.00 | 2 of 10 |
| BatchedQuickMC 1 | 10.52 | 1.54 | 2.84 | 11219 | 254 | 5.15 | 1.71 | 2 of 10 |
| default keys, fused | 6.78 | 2.39 | 0.76 | 39803 | 19 | 5.16 | 0.00 | 2 of 10 |
| BatchedQuickMC 1, fused | 7.88 | 2.05 | 0.21 | 11219 | 19 | 5.15 | 1.71 | 2 of 10 |
<!-- /table:s2_step0_search_e0 -->

<!-- table:s2_step0_search_three -->
| arm | wall | vs_default | evaluate | evaluate_n | us_per_eval | screen | bqmc | wrong |
|---|---|---|---|---|---|---|---|---|
| default keys | 5.44 | 1.00 | 4.63 | 8311 | 557 | 0.63 | 0.00 | 1 of 1 |
| BatchedQuickMC 1 | 4.53 | 1.20 | 0.99 | 1812 | 548 | 0.63 | 2.76 | 0 of 1 |
| default keys, fused | 1.70 | 3.20 | 0.92 | 8084 | 114 | 0.63 | 0.00 | 1 of 1 |
| BatchedQuickMC 1, fused | 3.67 | 1.48 | 0.20 | 1812 | 110 | 0.63 | 2.72 | 0 of 1 |
<!-- /table:s2_step0_search_three -->
