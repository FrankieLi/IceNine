<!-- table:level0_levels -->
| level | n_fz_cpp | cand_cpp | n_fz_py | cand_py | n_fz_probe | cand_probe |
|---|---|---|---|---|---|---|
| 0 | 4886 | 8151 | 4886 | 10190 | 4886 | 8110 |
| 1 | 2037 | 2525 | 2547 | 3504 | 2027 | 2487 |
| 2 | 631 | 691 | 876 | 1114 | 621 | 697 |
| 3 | 172 | 195 | 278 | 341 | 174 | 204 |
<!-- /table:level0_levels -->

<!-- table:level0_peak_patterns -->
| pattern | n_peaks | cpp_peak_overlap_sum |
|---|---|---|
| det0 bright=False, det1 lit=True, det1 bright=False | 29 | 0 |
| det0 bright=False, det1 lit=True, det1 bright=True | 2 | 0 |
| det0 bright=True, det1 lit=True, det1 bright=False | 6 | 0 |
| det0 bright=True, det1 lit=True, det1 bright=True | 1 | 1 |
| neither detector lit | 26 | 0 |
| not on a detector (no projection hit) | 12 | 0 |
<!-- /table:level0_peak_patterns -->

<!-- table:level0_s3_cpp_inputs -->
| voxel | n_cliques_with_ge2 | same_kept_set | same_kept_order |
|---|---|---|---|
| 15901 | 3537 | 3537 | 3537 |
<!-- /table:level0_s3_cpp_inputs -->

<!-- table:level0_settings -->
| level | n_fz_cpp | n_fz_py | n_local_cpp | n_local_py | nqmax_cpp | n_recip_cpp | n_recip_py | diam_deg_cpp | diam_deg_py | n_sym_cpp | n_sym_py |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 4886 | 9 | 9 | 5 | 26 | 26 | 5.000 | 5.000 | 24 | 24 |
| 1 | 2037 | 2547 | 9 | 9 | 6 | 50 | 50 | 3.333 | 3.333 | 24 | 24 |
| 2 | 631 | 876 | 9 | 9 | 7 | 64 | 64 | 2.222 | 2.222 | 24 | 24 |
| 3 | 172 | 278 | 9 | 9 | 8 | 112 | 112 | 1.481 | 1.481 | 24 | 24 |
<!-- /table:level0_settings -->

<!-- table:level0_stages -->
| python | screened_pairs_cpp | screened_pairs_py | fz_identical_screen | peak_overlap_equal | recost_median_diff | recost_within_1e6 | cliques_cpp | cliques_py | kept_cpp | kept_py |
|---|---|---|---|---|---|---|---|---|---|---|
| library default | 12326 | 35267 | 7 | 0.273 | -0.03 | 0.104 | 4491 | 4886 | 8151 | 10190 |
| C++ overlap rule | 12326 | 12485 | 3341 | 0.940 | 0 | 0.849 | 4491 | 4513 | 8151 | 8206 |
| C++ rule + centroid | 12326 | 12338 | 4466 | 0.986 | 0 | 0.853 | 4491 | 4480 | 8151 | 8110 |
<!-- /table:level0_stages -->

<!-- table:level0_voxels -->
| voxel | cpp | py_default | py_probe | ratio_default | ratio_probe |
|---|---|---|---|---|---|
| 15901 | 8,151 | 10,190 | 8,110 | 1.250 | 0.995 |
| 17034 | 8,284 | 10,382 | 8,283 | 1.253 | 1.000 |
| 18558 | 7,026 | 10,271 | 7,048 | 1.462 | 1.003 |
| 19889 | 7,673 | 10,308 | 7,683 | 1.343 | 1.001 |
| 22867 | 8,234 | 10,236 | 8,233 | 1.243 | 1.000 |
| 5242 | 6,674 | 10,213 | 6,674 | 1.530 | 1.000 |
<!-- /table:level0_voxels -->

<!-- table:level0_six_candidates -->
| voxel | l0_cpp | l0_before | l0_rule | all_cpp | all_before | all_rule |
|---|---|---|---|---|---|---|
| 15901 | 8,151 | 10,190 | 8,110 | 11,562 | 15,149 | 11,498 |
| 17034 | 8,284 | 10,382 | 8,283 | 11,609 | 15,483 | 11,689 |
| 18558 | 7,026 | 10,271 | 7,048 | 9,973 | 15,349 | 10,004 |
| 19889 | 7,673 | 10,308 | 7,683 | 10,820 | 15,341 | 10,877 |
| 22867 | 8,234 | 10,236 | 8,233 | 11,549 | 15,277 | 11,613 |
| 5242 | 6,674 | 10,213 | 6,674 | 9,461 | 15,144 | 9,452 |
<!-- /table:level0_six_candidates -->

<!-- table:level0_six_evals -->
| voxel | evals_cpp | evals_before | evals_rule | evals_ratio_before | evals_ratio_rule | wall_cpp_s | wall_before_s | wall_rule_s | label_rule | err_before | err_rule |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 251,188 | 328,519 | 237,739 | 1.31 | 0.95 | 21.2 | 176 | 125 | single-process | 0.0170 | 0.0291 |
| 17034 | 251,673 | 339,561 | 246,127 | 1.35 | 0.98 | 21.0 | 179 | 125 | single-process, contended | 0.0284 | 0.0639 |
| 18558 | 219,790 | 329,367 | 209,728 | 1.50 | 0.95 | 16.2 | 171 | 106 | single-process | 0.0048 | 0.0034 |
| 19889 | 236,372 | 330,209 | 226,453 | 1.40 | 0.96 | 18.4 | 176 | 115 | single-process | 0.0210 | 0.0212 |
| 22867 | 250,852 | 335,687 | 244,186 | 1.34 | 0.97 | 19.1 | 174 | 125 | single-process | 0.0263 | 0.0530 |
| 5242 | 209,927 | 330,577 | 204,646 | 1.57 | 0.97 | 10.5 | 171 | 102 | single-process | 0.0493 | 0.0281 |
<!-- /table:level0_six_evals -->

<!-- table:cpprule_e0_on -->
| variant | n | wrong_before | wrong_after | broken | fixed | both_wrong | med_both_before | med_both_after | evals_before | evals_after | basin_miss_before | basin_miss_after |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40/200 [0.150, 0.261] | 61/200 [0.245, 0.372] | 38 | 17 | 23 | 0.0102 | 0.0187 | 48442 | 47174 | 4 | 6 |
| realistic | 200 | 25/200 [0.086, 0.178] | 42/200 [0.159, 0.272] | 25 | 8 | 17 | 0.0079 | 0.0157 | 51058 | 48731 | 2 | 10 |
<!-- /table:cpprule_e0_on -->

<!-- table:cpprule_e0_rule -->
| variant | n | wrong_before | wrong_after | broken | fixed | both_wrong | med_both_before | med_both_after | evals_before | evals_after | basin_miss_before | basin_miss_after |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40/200 [0.150, 0.261] | 63/200 [0.255, 0.382] | 39 | 16 | 24 | 0.0093 | 0.0131 | 48442 | 47067 | 4 | 8 |
| realistic | 200 | 25/200 [0.086, 0.178] | 54/200 [0.213, 0.335] | 36 | 7 | 18 | 0.0080 | 0.0118 | 51058 | 48487 | 2 | 9 |
<!-- /table:cpprule_e0_rule -->

<!-- table:e0_cpp_vs_python -->
| variant | n | py_off | py_on | cpp | newly_wrong | cpp_wrong_among_newly | cpp_only_vs_off | off_only_vs_cpp | cpp_only_vs_on | on_only_vs_cpp |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40 | 61 | 71 [0.292, 0.423] | 38 | 24 | 47 | 16 | 28 | 18 |
| realistic | 200 | 25 | 42 | 52 [0.204, 0.325] | 25 | 11 | 38 | 11 | 28 | 18 |
<!-- /table:e0_cpp_vs_python -->

<!-- table:e0_fairness_broken -->
| setting | variant | n | off_wrong | on_wrong | broken | fixed | both_wrong | cpp_s0_wrong_on_broken | cpp_s1_wrong_on_broken | cpp_s2_wrong_on_broken |
|---|---|---|---|---|---|---|---|---|---|---|
| min_sin_eta 0.3 | clean | 200 | 40 | 61 | 38 | 17 | 23 | 24 | 13 | 17 |
| min_sin_eta 0.3 | realistic | 200 | 25 | 42 | 25 | 8 | 17 | 11 | 10 | 16 |
| min_sin_eta 0 | clean | 200 | 39 | 56 | 34 | 17 | 22 | 18 | 13 | 19 |
| min_sin_eta 0 | realistic | 200 | 18 | 56 | 44 | 6 | 12 | 22 | 19 | 22 |
<!-- /table:e0_fairness_broken -->

<!-- table:e0_fairness_counts -->
| variant | n | py_off | py_on | py_off_ms0 | py_on_ms0 | cpp_s0 | cpp_s1 | cpp_s2 | cpp_mean | cpp_wrong_all_seeds | cpp_wrong_any_seed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40 | 61 | 39 | 56 | 71 | 56 | 66 | 64.3 | 28 | 105 |
| realistic | 200 | 25 | 42 | 18 | 56 | 52 | 46 | 52 | 50.0 | 25 | 83 |
<!-- /table:e0_fairness_counts -->

<!-- table:e0_debug_classes -->
| cls | py_off | py_on | cpp |
|---|---|---|---|
| never | 0 | 4 | 4 |
| pruned 0 | 0 | 25 | 15 |
| pruned 1 | 0 | 23 | 11 |
| pruned 2 | 0 | 11 | 5 |
| right | 63 | 0 | 28 |
<!-- /table:e0_debug_classes -->

<!-- table:e0_debug_level0_match -->
| py | n_py_median | n_cpp_median | cpp_cands_matched_frac_median | median_abs_cost_diff |
|---|---|---|---|---|
| on | 125 | 124 | 0.978 | 0.00e+00 |
| off | 178 | 124 | 0.952 | 0.00e+00 |
<!-- /table:e0_debug_level0_match -->

<!-- table:e0_debug_levels -->
| run | level | n | n_cand | closest_deg | rank | n_keep | survives | cost_closest | cost_best | within3 |
|---|---|---|---|---|---|---|---|---|---|---|
| py_off | 0 | 63 | 178 | 1.12 | 9 | 44 | 53 | 0.962 | 0.927 | 63 |
| py_off | 1 | 63 | 46 | 0.38 | 0 | 11 | 63 | 0.833 | 0.827 | 63 |
| py_off | 2 | 63 | 12 | 0.20 | 0 | 3 | 63 | 0.679 | 0.665 | 63 |
| py_off | 3 | 63 | 3 | 0.16 | 0 | 3 | 63 | 0.526 | 0.526 | 63 |
| py_on | 0 | 63 | 125 | 1.56 | 29 | 31 | 28 | 0.981 | 0.937 | 59 |
| py_on | 1 | 63 | 31 | 2.26 | 15 | 7 | 17 | 0.967 | 0.907 | 35 |
| py_on | 2 | 63 | 7 | 7.22 | 4 | 1 | 2 | 0.951 | 0.861 | 12 |
| py_on | 3 | 63 | 1 | 38.19 | 0 | 1 | 63 | 0.875 | 0.843 | 0 |
| cpp | 0 | 63 | 124 | 1.27 | 18 | 31 | 41 | 0.976 | 0.935 | 59 |
| cpp | 1 | 63 | 32 | 1.09 | 5 | 8 | 37 | 0.954 | 0.896 | 44 |
| cpp | 2 | 63 | 8 | 2.21 | 2 | 2 | 28 | 0.940 | 0.814 | 33 |
| cpp | 3 | 63 | 2 | 18.94 | 0 | 2 | 63 | 0.797 | 0.782 | 28 |
<!-- /table:e0_debug_levels -->

<!-- table:level0_cost_offset -->
| rule | n_within_1e6 | n_within_1e3 | mean_signed | median_signed | n_python_lower | n_python_higher | corr |
|---|---|---|---|---|---|---|---|
| off (library default) | 8 | 10 | -0.01311 | -1.24e-02 | 151 | 1 | 0.99896 |
| on | 143 | 154 | 0.00005 | -3.95e-10 | 2 | 15 | 0.99997 |
<!-- /table:level0_cost_offset -->

<!-- table:e0_cost_check -->
| py_min_sin_eta | py_rule | n | median_abs_dq | max_abs_dq | median_abs_dconf | max_abs_dconf | same_peak_counts |
|---|---|---|---|---|---|---|---|
| 0.3 | off | 40 | 0.0093 | 0.0391 | 0.0892 | 0.2405 | 1 |
| 0.3 | on | 40 | 0.0029 | 0.0198 | 0.0058 | 0.0212 | 3 |
| 0.0 | off | 40 | 0.0138 | 0.0446 | 0.0830 | 0.2273 | 7 |
| 0.0 | on | 40 | 0.0000 | 0.0161 | 0.0000 | 0.0161 | 29 |
<!-- /table:e0_cost_check -->
