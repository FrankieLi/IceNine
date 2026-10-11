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
| variant | n | wrong_before | wrong_after | only_after | only_before | p | med_both_before | med_both_after | evals_before | evals_after | basin_miss_before | basin_miss_after |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40/200 [0.150, 0.261] | 61/200 [0.245, 0.372] | 38 | 17 | 0.0065 | 0.0102 | 0.0187 | 48442 | 47174 | 4 | 6 |
| realistic | 200 | 25/200 [0.086, 0.178] | 42/200 [0.159, 0.272] | 25 | 8 | 0.0046 | 0.0079 | 0.0157 | 51058 | 48731 | 2 | 10 |
<!-- /table:cpprule_e0_on -->

<!-- table:cpprule_e0_rule -->
| variant | n | wrong_before | wrong_after | only_after | only_before | p | med_both_before | med_both_after | evals_before | evals_after | basin_miss_before | basin_miss_after |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40/200 [0.150, 0.261] | 63/200 [0.255, 0.382] | 39 | 16 | 0.0027 | 0.0093 | 0.0131 | 48442 | 47067 | 4 | 8 |
| realistic | 200 | 25/200 [0.086, 0.178] | 54/200 [0.213, 0.335] | 36 | 7 | 9e-06 | 0.0080 | 0.0118 | 51058 | 48487 | 2 | 9 |
<!-- /table:cpprule_e0_rule -->

<!-- table:level0_cost_offset -->
| rule | n_within_1e6 | n_within_1e3 | mean_signed | median_signed | n_python_lower | n_python_higher | corr |
|---|---|---|---|---|---|---|---|
| off (library default) | 8 | 10 | -0.01311 | -1.24e-02 | 151 | 1 | 0.99896 |
| on | 143 | 154 | 0.00005 | -3.95e-10 | 2 | 15 | 0.99997 |
<!-- /table:level0_cost_offset -->
