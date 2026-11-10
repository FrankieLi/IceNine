<!-- table:evalA_cmp_levels_cpp_edge_interp_batched_v15901 -->
| level | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | kept_py | kept_cpp |
|---|---|---|---|---|---|---|---|---|
| 0 | 4,886 | 4,886 | 12,325 | 12,326 | 8,149 | 8,151 | 2,037 | 2,037 |
| 1 | 2,037 | 2,037 | 11,035 | 10,852 | 2,465 | 2,525 | 616 | 631 |
| 2 | 616 | 631 | 4,131 | 4,207 | 695 | 691 | 173 | 172 |
| 3 | 173 | 172 | 1,468 | 1,488 | 198 | 195 | 49 | 48 |
<!-- /table:evalA_cmp_levels_cpp_edge_interp_batched_v15901 -->

<!-- table:evalA_cmp_levels_cpp_edge_interp_v15901 -->
| level | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | kept_py | kept_cpp |
|---|---|---|---|---|---|---|---|---|
| 0 | 4,886 | 4,886 | 12,325 | 12,326 | 8,149 | 8,151 | 2,037 | 2,037 |
| 1 | 2,037 | 2,037 | 11,035 | 10,852 | 2,465 | 2,525 | 616 | 631 |
| 2 | 616 | 631 | 4,131 | 4,207 | 695 | 691 | 173 | 172 |
| 3 | 173 | 172 | 1,468 | 1,488 | 198 | 195 | 49 | 48 |
<!-- /table:evalA_cmp_levels_cpp_edge_interp_v15901 -->

<!-- table:evalA_cmp_levels_cpp_edge_v15901 -->
| level | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | kept_py | kept_cpp |
|---|---|---|---|---|---|---|---|---|
| 0 | 4,886 | 4,886 | 12,298 | 12,326 | 8,107 | 8,151 | 2,026 | 2,037 |
| 1 | 2,026 | 2,037 | 10,920 | 10,852 | 2,481 | 2,525 | 620 | 631 |
| 2 | 620 | 631 | 4,187 | 4,207 | 683 | 691 | 170 | 172 |
| 3 | 170 | 172 | 1,466 | 1,488 | 191 | 195 | 47 | 48 |
<!-- /table:evalA_cmp_levels_cpp_edge_v15901 -->

<!-- table:evalA_cmp_levels_cpp_v15901 -->
| level | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | kept_py | kept_cpp |
|---|---|---|---|---|---|---|---|---|
| 0 | 4,886 | 4,886 | 12,338 | 12,326 | 8,110 | 8,151 | 2,027 | 2,037 |
| 1 | 2,027 | 2,037 | 10,909 | 10,852 | 2,487 | 2,525 | 621 | 631 |
| 2 | 621 | 631 | 4,181 | 4,207 | 697 | 691 | 174 | 172 |
| 3 | 174 | 172 | 1,492 | 1,488 | 204 | 195 | 51 | 48 |
<!-- /table:evalA_cmp_levels_cpp_v15901 -->

<!-- table:evalA_e0_interp_paired -->
| new | base | n | wrong_base | wrong_new | fixed | broken | both_wrong | evals_base | evals_new | l0_base | l0_new | basin_miss_base | basin_miss_new | runtime_base | runtime_new |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| edge | on | 200 | 56/200 [0.222, 0.346] | 62/200 [0.250, 0.377] | 21 | 27 | 35 | 47,230 | 47,216 | 125 | 125 | 6 | 6 | 24.2 | 24.8 |
| interp | on | 200 | 56/200 [0.222, 0.346] | 65/200 [0.264, 0.393] | 20 | 29 | 36 | 47,230 | 47,179 | 125 | 122 | 6 | 7 | 24.2 | 3.4 |
| interp | edge | 200 | 62/200 [0.250, 0.377] | 65/200 [0.264, 0.393] | 20 | 23 | 42 | 47,216 | 47,179 | 125 | 122 | 6 | 7 | 24.8 | 3.4 |
| lenint | off | 200 | 41/200 [0.155, 0.266] | 40/200 [0.150, 0.261] | 21 | 20 | 20 | 48,700 | 48,633 | 182 | 180 | 4 | 4 | 26.0 | 4.4 |
<!-- /table:evalA_e0_interp_paired -->

<!-- table:evalA_interp_by_j -->
| j | local_rotation_deg | n | same_pass_pct | same_counts_pct | omega_err_median_of_max_deg | omega_err_worst_deg | pct_with_exact_peak_missing | pct_with_extra_peak |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.000 | 102 | 100.0 | 100.0 | 0.0000 | 0.0000 | 0.0 | 0.0 |
| 1 | 2.837 | 92 | 96.7 | 93.5 | 0.0896 | 6.5947 | 2.2 | 1.1 |
| 2 | 2.837 | 98 | 100.0 | 98.0 | 0.0813 | 174.5162 | 3.1 | 2.0 |
| 3 | 2.837 | 109 | 99.1 | 98.2 | 0.1278 | 169.9071 | 3.7 | 2.8 |
| 4 | 2.837 | 94 | 96.8 | 93.6 | 0.1174 | 170.1070 | 2.1 | 2.1 |
| 5 | 2.837 | 134 | 99.3 | 96.3 | 0.0945 | 174.5162 | 2.2 | 2.2 |
| 6 | 2.837 | 104 | 100.0 | 95.2 | 0.1089 | 133.9467 | 4.8 | 3.8 |
| 7 | 2.837 | 73 | 97.3 | 97.3 | 0.1066 | 18.5853 | 1.4 | 0.0 |
| 8 | 2.837 | 94 | 100.0 | 94.7 | 0.1040 | 170.1794 | 6.4 | 3.2 |
<!-- /table:evalA_interp_by_j -->

<!-- table:evalA_interp_by_theta -->
| theta_deg | n | same_pass_pct | same_counts_pct | same_counts_lo | same_counts_hi | omega_err_median_of_max_deg | omega_err_p90_of_max_deg | omega_err_worst_deg | pct_with_exact_peak_missing | pct_with_extra_peak |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0000 | 0.0005 | 0.602 | 0.0 | 0.0 |
| 0.10 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0001 | 0.0011 | 0.631 | 0.0 | 0.0 |
| 0.20 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0004 | 0.0033 | 0.074 | 0.0 | 0.0 |
| 0.50 | 75 | 100.0 | 98.7 | 92.8 | 99.8 | 0.0021 | 0.0818 | 174.464 | 1.3 | 1.3 |
| 1.00 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0113 | 0.1912 | 7.303 | 0.0 | 0.0 |
| 2.00 | 75 | 100.0 | 98.7 | 92.8 | 99.8 | 0.0428 | 1.3502 | 170.192 | 1.3 | 1.3 |
| 3.00 | 75 | 98.7 | 96.0 | 88.9 | 98.6 | 0.1182 | 0.4056 | 173.865 | 2.7 | 1.3 |
| 5.00 | 75 | 100.0 | 97.3 | 90.8 | 99.3 | 0.3495 | 1.6521 | 8.550 | 4.0 | 4.0 |
<!-- /table:evalA_interp_by_theta -->

<!-- table:evalA_interp_vs_cpp -->
| group | n | py_interp_pass_matches_cpp | py_interp_counts_match_cpp | py_exact_pass_matches_cpp_interp |
|---|---|---|---|---|
| pair_random | 900 | 900 [0.9957, 1.0000] | 900 | 890 |
| sweep | 600 | 600 [0.9936, 1.0000] | 600 | 599 |
| pair_differ_pass | 442 | 442 [0.9914, 1.0000] | 442 | 40 |
| pair_differ_count | 58 | 58 [0.9379, 1.0000] | 56 | 58 |
| ALL | 2000 | 2000 [0.9981, 1.0000] | 1998 | 1587 |
<!-- /table:evalA_interp_vs_cpp -->

<!-- table:evalA_ledger_mean -->
| arm | n | screen | recost | quick | find | variance | final | total | wall_s | us_per_eval |
|---|---|---|---|---|---|---|---|---|---|---|
| base | 6 | 77,510 | 67,084 | 183,148 | 2,885 | 1,418 | 1 | 332,046 | 173 | 491 |
| cpp | 6 | 67,978 | 26,148 | 130,266 | 2,354 | 1,398 | 1 | 228,146 | 116 | 486 |
| cpp_edge | 6 | 67,968 | 26,077 | 130,184 | 198 | 1,455 | 1 | 225,883 | 119 | 506 |
| cpp_edge_batched | 6 | 67,968 | 26,077 | 130,184 | 198 | 1,455 | 1 | 225,883 | 93 | 388 |
| cpp_edge_interp | 6 | 67,911 | 26,016 | 129,894 | 182 | 1,644 | 1 | 225,648 | 97 | 410 |
| cpp_edge_interp_batched | 6 | 67,911 | 26,016 | 129,894 | 182 | 1,644 | 1 | 225,648 | 92 | 386 |
| lenient_interp_batched | 6 | 77,540 | 67,050 | 183,304 | 238 | 1,462 | 1 | 329,595 | 150 | 426 |
| r2 | 6 | 69,416 | 30,522 | 138,052 | 2,444 | 1,398 | 1 | 241,832 | 132 | 525 |
<!-- /table:evalA_ledger_mean -->

<!-- table:evalA_ledger_stages -->
| arm | voxel | screen | recost | quick | find | variance | final | total | wall_s | label | err |
|---|---|---|---|---|---|---|---|---|---|---|---|
| base | 15901 | 77,139 | 67,582 | 181,068 | 530 | 1,299 | 1 | 327,619 | 167 | single-worker | 0.038 |
| base | 17034 | 78,066 | 68,860 | 186,312 | 5,114 | 2,069 | 1 | 340,422 | 179 | single-worker, contended | 0.022 |
| base | 18558 | 77,580 | 65,831 | 183,720 | 1,253 | 1,530 | 1 | 329,915 | 172 | single-worker, contended | 0.005 |
| base | 19889 | 77,643 | 67,062 | 183,744 | 129 | 1,145 | 1 | 329,724 | 172 | single-worker, contended | 0.015 |
| base | 22867 | 77,220 | 67,930 | 181,452 | 5,144 | 1,299 | 1 | 333,046 | 173 | single-worker, contended | 0.012 |
| base | 5242 | 77,409 | 65,241 | 182,592 | 5,142 | 1,167 | 1 | 331,552 | 173 | single-worker | 0.060 |
| cpp | 15901 | 69,372 | 28,920 | 137,976 | 193 | 1,277 | 1 | 237,739 | 121 | single-worker | 0.029 |
| cpp | 17034 | 69,849 | 29,641 | 140,268 | 4,376 | 1,992 | 1 | 246,127 | 125 | single-worker, contended | 0.064 |
| cpp | 18558 | 66,096 | 21,922 | 120,048 | 208 | 1,453 | 1 | 209,728 | 106 | single-worker, contended | 0.003 |
| cpp | 19889 | 68,040 | 26,626 | 130,524 | 161 | 1,101 | 1 | 226,453 | 116 | single-worker, contended | 0.021 |
| cpp | 22867 | 69,651 | 29,663 | 139,356 | 4,205 | 1,310 | 1 | 244,186 | 124 | single-worker, contended | 0.053 |
| cpp | 5242 | 64,863 | 20,119 | 113,424 | 4,984 | 1,255 | 1 | 204,646 | 103 | single-worker, contended | 0.028 |
| cpp_edge | 15901 | 69,318 | 28,871 | 137,544 | 208 | 1,948 | 1 | 237,890 | 131 | single-worker, contended | 0.016 |
| cpp_edge | 17034 | 69,759 | 29,632 | 139,752 | 161 | 1,321 | 1 | 240,626 | 127 | single-worker, contended | 0.044 |
| cpp_edge | 18558 | 66,159 | 21,948 | 120,468 | 401 | 1,266 | 1 | 210,243 | 110 | single-worker, contended | 0.005 |
| cpp_edge | 19889 | 68,004 | 26,527 | 130,284 | 161 | 1,398 | 1 | 226,375 | 118 | single-worker, contended | 0.034 |
| cpp_edge | 22867 | 69,660 | 29,504 | 139,428 | 161 | 2,025 | 1 | 240,779 | 126 | single-worker | 0.038 |
| cpp_edge | 5242 | 64,908 | 19,981 | 113,628 | 97 | 771 | 1 | 199,386 | 104 | single-worker, contended | 0.012 |
| cpp_edge_batched | 15901 | 69,318 | 28,871 | 137,544 | 208 | 1,948 | 1 | 237,890 | 98 | single-worker, contended | 0.016 |
| cpp_edge_batched | 17034 | 69,759 | 29,632 | 139,752 | 161 | 1,321 | 1 | 240,626 | 103 | single-worker, contended | 0.044 |
| cpp_edge_batched | 18558 | 66,159 | 21,948 | 120,468 | 401 | 1,266 | 1 | 210,243 | 87 | single-worker | 0.005 |
| cpp_edge_batched | 19889 | 68,004 | 26,527 | 130,284 | 161 | 1,398 | 1 | 226,375 | 91 | single-worker, contended | 0.034 |
| cpp_edge_batched | 22867 | 69,660 | 29,504 | 139,428 | 161 | 2,025 | 1 | 240,779 | 99 | single-worker, contended | 0.038 |
| cpp_edge_batched | 5242 | 64,908 | 19,981 | 113,628 | 97 | 771 | 1 | 199,386 | 78 | single-worker | 0.012 |
| cpp_edge_interp | 15901 | 69,408 | 28,959 | 138,084 | 193 | 2,223 | 1 | 238,868 | 103 | single-worker, contended | 0.013 |
| cpp_edge_interp | 17034 | 69,714 | 29,476 | 139,560 | 208 | 2,245 | 1 | 241,204 | 109 | single-worker, contended | 0.014 |
| cpp_edge_interp | 18558 | 66,177 | 21,849 | 120,480 | 208 | 1,255 | 1 | 209,970 | 88 | single-worker | 0.005 |
| cpp_edge_interp | 19889 | 67,878 | 26,347 | 129,648 | 193 | 1,002 | 1 | 225,069 | 99 | single-worker, contended | 0.026 |
| cpp_edge_interp | 22867 | 69,453 | 29,438 | 138,360 | 161 | 2,168 | 1 | 239,581 | 102 | single-worker | 0.022 |
| cpp_edge_interp | 5242 | 64,836 | 20,028 | 113,232 | 129 | 969 | 1 | 199,195 | 85 | single-worker | 0.056 |
| cpp_edge_interp_batched | 15901 | 69,408 | 28,959 | 138,084 | 193 | 2,223 | 1 | 238,868 | 98 | single-worker, contended | 0.013 |
| cpp_edge_interp_batched | 17034 | 69,714 | 29,476 | 139,560 | 208 | 2,245 | 1 | 241,204 | 99 | single-worker | 0.014 |
| cpp_edge_interp_batched | 18558 | 66,177 | 21,849 | 120,480 | 208 | 1,255 | 1 | 209,970 | 84 | single-worker | 0.005 |
| cpp_edge_interp_batched | 19889 | 67,878 | 26,347 | 129,648 | 193 | 1,002 | 1 | 225,069 | 94 | single-worker | 0.026 |
| cpp_edge_interp_batched | 22867 | 69,453 | 29,438 | 138,360 | 161 | 2,168 | 1 | 239,581 | 99 | single-worker | 0.022 |
| cpp_edge_interp_batched | 5242 | 64,836 | 20,028 | 113,232 | 129 | 969 | 1 | 199,195 | 78 | single-worker | 0.056 |
| lenient_interp_batched | 15901 | 77,256 | 67,586 | 181,608 | 305 | 1,959 | 1 | 328,715 | 152 | single-worker | 0.015 |
| lenient_interp_batched | 17034 | 77,841 | 68,523 | 184,932 | 208 | 1,629 | 1 | 333,134 | 149 | single-worker | 0.024 |
| lenient_interp_batched | 18558 | 77,445 | 65,605 | 182,772 | 466 | 1,332 | 1 | 327,621 | 149 | single-worker | 0.004 |
| lenient_interp_batched | 19889 | 77,697 | 67,071 | 184,092 | 161 | 1,387 | 1 | 330,409 | 147 | single-worker | 0.030 |
| lenient_interp_batched | 22867 | 77,472 | 68,191 | 183,120 | 97 | 1,420 | 1 | 330,301 | 152 | single-worker | 0.042 |
| lenient_interp_batched | 5242 | 77,526 | 65,324 | 183,300 | 193 | 1,046 | 1 | 327,390 | 151 | single-worker | 0.022 |
| r2 | 15901 | 70,587 | 33,610 | 144,336 | 387 | 1,574 | 1 | 250,495 | 135 | single-worker, contended | 0.014 |
| r2 | 17034 | 70,794 | 34,155 | 145,416 | 4,732 | 1,673 | 1 | 256,771 | 144 | single-worker, contended | 0.025 |
| r2 | 18558 | 67,725 | 25,638 | 128,964 | 337 | 1,079 | 1 | 223,744 | 125 | single-worker | 0.005 |
| r2 | 19889 | 69,939 | 31,476 | 140,844 | 161 | 1,156 | 1 | 243,577 | 133 | single-worker | 0.033 |
| r2 | 22867 | 70,938 | 34,334 | 146,400 | 4,318 | 1,838 | 1 | 257,829 | 140 | single-worker | 0.061 |
| r2 | 5242 | 66,510 | 23,918 | 122,352 | 4,728 | 1,068 | 1 | 218,577 | 116 | single-worker | 0.023 |
<!-- /table:evalA_ledger_stages -->

<!-- table:evalA_levels_base_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 35232 | 35232 | 10187 | 122244 | 2546 | 7641 | 17.0 | 17.8 | 37.8 | 65.1 |
| 1 | 2546 | 9 | 22914 | 22180 | 22180 | 3461 | 41532 | 865 | 2596 | 9.7 | 11.1 | 22.0 | 21.8 |
| 2 | 865 | 9 | 7785 | 7705 | 7705 | 1098 | 13176 | 274 | 824 | 3.5 | 3.8 | 7.7 | 6.9 |
| 3 | 274 | 9 | 2466 | 2465 | 2465 | 343 | 4116 | 85 | 258 | 1.3 | 1.2 | 2.6 | 2.1 |
<!-- /table:evalA_levels_base_v15901 -->

<!-- table:evalA_levels_cpp_edge_batched_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 12298 | 12298 | 8107 | 97284 | 2026 | 6081 | 1.1 | 6.6 | 8.5 | 54.8 |
| 1 | 2026 | 9 | 18234 | 10920 | 10920 | 2481 | 29772 | 620 | 1861 | 0.7 | 5.8 | 7.1 | 16.9 |
| 2 | 620 | 9 | 5580 | 4187 | 4187 | 683 | 8196 | 170 | 513 | 0.3 | 2.2 | 2.7 | 4.6 |
| 3 | 170 | 9 | 1530 | 1466 | 1466 | 191 | 2292 | 47 | 144 | 0.1 | 0.8 | 1.0 | 1.3 |
<!-- /table:evalA_levels_cpp_edge_batched_v15901 -->

<!-- table:evalA_levels_cpp_edge_interp_batched_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 12325 | 12325 | 8149 | 97788 | 2037 | 6112 | 1.1 | 6.6 | 8.4 | 55.1 |
| 1 | 2037 | 9 | 18333 | 11035 | 11035 | 2465 | 29580 | 616 | 1849 | 0.8 | 5.9 | 7.2 | 16.6 |
| 2 | 616 | 9 | 5544 | 4131 | 4131 | 695 | 8340 | 173 | 522 | 0.3 | 2.2 | 2.7 | 4.7 |
| 3 | 173 | 9 | 1557 | 1468 | 1468 | 198 | 2376 | 49 | 149 | 0.1 | 0.8 | 1.0 | 1.3 |
<!-- /table:evalA_levels_cpp_edge_interp_batched_v15901 -->

<!-- table:evalA_levels_cpp_edge_interp_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 12325 | 12325 | 8149 | 97788 | 2037 | 6112 | 3.8 | 6.6 | 11.2 | 54.9 |
| 1 | 2037 | 9 | 18333 | 11035 | 11035 | 2465 | 29580 | 616 | 1849 | 2.0 | 5.9 | 8.4 | 16.5 |
| 2 | 616 | 9 | 5544 | 4131 | 4131 | 695 | 8340 | 173 | 522 | 0.7 | 2.2 | 3.0 | 4.7 |
| 3 | 173 | 9 | 1557 | 1468 | 1468 | 198 | 2376 | 49 | 149 | 0.3 | 0.9 | 1.2 | 1.4 |
<!-- /table:evalA_levels_cpp_edge_interp_v15901 -->

<!-- table:evalA_levels_cpp_edge_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 12298 | 12298 | 8107 | 97284 | 2026 | 6081 | 19.2 | 6.9 | 27.0 | 57.1 |
| 1 | 2026 | 9 | 18234 | 10920 | 10920 | 2481 | 29772 | 620 | 1861 | 8.7 | 6.1 | 15.3 | 17.3 |
| 2 | 620 | 9 | 5580 | 4187 | 4187 | 683 | 8196 | 170 | 513 | 2.8 | 2.3 | 5.3 | 4.8 |
| 3 | 170 | 9 | 1530 | 1466 | 1466 | 191 | 2292 | 47 | 144 | 0.9 | 0.8 | 1.7 | 1.3 |
<!-- /table:evalA_levels_cpp_edge_v15901 -->

<!-- table:evalA_levels_cpp_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 12338 | 12338 | 8110 | 97320 | 2027 | 6083 | 17.4 | 6.4 | 24.6 | 53.1 |
| 1 | 2027 | 9 | 18243 | 10909 | 10909 | 2487 | 29844 | 621 | 1866 | 7.9 | 5.6 | 14.1 | 16.4 |
| 2 | 621 | 9 | 5589 | 4181 | 4181 | 697 | 8364 | 174 | 523 | 2.5 | 2.1 | 4.9 | 4.5 |
| 3 | 174 | 9 | 1566 | 1492 | 1492 | 204 | 2448 | 51 | 153 | 0.8 | 0.8 | 1.7 | 1.3 |
<!-- /table:evalA_levels_cpp_v15901 -->

<!-- table:evalA_levels_lenient_interp_batched_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 35186 | 35186 | 10206 | 122472 | 2551 | 7655 | 1.1 | 19.6 | 23.7 | 71.6 |
| 1 | 2551 | 9 | 22959 | 22157 | 22157 | 3475 | 41700 | 868 | 2607 | 1.0 | 12.4 | 14.7 | 24.1 |
| 2 | 868 | 9 | 7812 | 7732 | 7732 | 1116 | 13392 | 279 | 837 | 0.4 | 4.2 | 5.1 | 7.7 |
| 3 | 279 | 9 | 2511 | 2511 | 2511 | 337 | 4044 | 84 | 253 | 0.2 | 1.4 | 1.7 | 2.3 |
<!-- /table:evalA_levels_lenient_interp_batched_v15901 -->

<!-- table:evalA_levels_r2_v15901 -->
| level | fz_in | local | screen | passed | recost | cands | quick | kept | dropped | screen_s | recost_s | discrete_s | quick_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4886 | 9 | 43974 | 14526 | 14526 | 8485 | 101820 | 2121 | 6364 | 18.5 | 7.9 | 27.4 | 58.4 |
| 1 | 2121 | 9 | 19089 | 12647 | 12647 | 2626 | 31512 | 656 | 1970 | 8.8 | 6.9 | 16.4 | 18.0 |
| 2 | 656 | 9 | 5904 | 4891 | 4891 | 722 | 8664 | 180 | 542 | 2.9 | 2.7 | 5.8 | 4.9 |
| 3 | 180 | 9 | 1620 | 1546 | 1546 | 195 | 2340 | 48 | 147 | 0.8 | 0.8 | 1.7 | 1.3 |
<!-- /table:evalA_levels_r2_v15901 -->

<!-- table:evalA_memory_batched_screen -->
| voxel | batched_screen | screen_chunk | rss_setup_mib | rss_peak_mib | evals | err_deg | wall_s | label |
|---|---|---|---|---|---|---|---|---|
| 15901 | off | 1 | 285.7 | 872.3 | 327,619 | 0.0378 | 173.0 | single-worker, contended |
| 15901 | on | 32 | 286.9 | 880.5 | 327,619 | 0.0378 | 145.4 | single-worker, contended |
<!-- /table:evalA_memory_batched_screen -->

<!-- table:evalA_parity_by_group -->
| group | n | EXS | EXS_edge | EXF | EXF_edge |
|---|---|---|---|---|---|
| truth | 3 | 2 | 3 | 3 | 3 |
| off_0.1 | 9 | 6 | 9 | 5 | 9 |
| off_1.0 | 9 | 6 | 9 | 3 | 9 |
| off_3.0 | 9 | 9 | 9 | 4 | 9 |
| random | 200 | 178 | 200 | 128 | 200 |
| pair_random | 900 | 815 | 900 | 553 | 900 |
| sweep | 600 | 547 | 600 | 382 | 600 |
| pair_differ_pass | 442 | 356 | 441 | 261 | 442 |
| pair_differ_count | 58 | 46 | 58 | 32 | 58 |
<!-- /table:evalA_parity_by_group -->

<!-- table:evalA_parity_causes -->
| arm | peaks | frame_below_range | frame_boundary | pixel | counts | omega_only |
|---|---|---|---|---|---|---|
| screen, Python default | 108196 | 283 | 1 | 0 | 0 | 32487 |
| screen, Python with CppEdgeRules | 108196 | 0 | 1 | 0 | 0 | 32586 |
| complete cost, Python default | 463012 | 856 | 1 | 269 | 311 | 137918 |
| complete cost, Python with CppEdgeRules | 463012 | 0 | 1 | 259 | 0 | 138368 |
<!-- /table:evalA_parity_causes -->

<!-- table:evalA_parity_exact -->
| arm | n | cost | counts | frames | qual |
|---|---|---|---|---|---|
| screen, Python default | 2230 | 1965 (88.1%) | 1965 (88.1%) | 1846 (82.8%) | n/a |
| screen, Python with CppEdgeRules | 2230 | 2229 (99.96%) | 2229 (99.96%) | 2187 (98.1%) | n/a |
| complete cost, Python default | 2230 | 1371 (61.5%) | 1124 (50.4%) | 887 (39.8%) | 2020 (90.6%) |
| complete cost, Python with CppEdgeRules | 2230 | 2230 (100.0%) | 1985 (89.0%) | 1690 (75.8%) | 2228 (99.9%) |
<!-- /table:evalA_parity_exact -->

<!-- table:evalA_pass_patterns -->
| group | pattern | n |
|---|---|---|
| pair_random | ---- | 664 |
| pair_random | PPPP | 226 |
| pair_random | P--- | 6 |
| pair_random | -PPP | 4 |
| pair_differ_pass | P--- | 214 |
| pair_differ_pass | -PPP | 187 |
| pair_differ_pass | --P- | 40 |
| pair_differ_pass | PP-- | 1 |
| pair_differ_count | PPPP | 58 |
<!-- /table:evalA_pass_patterns -->

<!-- table:evalA_passsets_l0 -->
| arm | edge_rules | pass_pairs | cpp_pass_pairs | fz_identical | pairs_pass_differs | overlap_equal | overlap_equal_pct |
|---|---|---|---|---|---|---|---|
| exact, one orientation per call (today) | off | 12,338 | 12,326 | 4,466 | 442 | 43,359 | 98.601 |
| exact, 32 FZ points per call (BatchedScreen) | off | 12,338 | 12,326 | 4,466 | 442 | 43,359 | 98.601 |
| interpolated, one FZ point per call (CppOmegaInterp) | off | 12,364 | 12,326 | 4,849 | 40 | 43,913 | 99.861 |
| interpolated, 32 FZ points per call (both) | off | 12,364 | 12,326 | 4,849 | 40 | 43,912 | 99.859 |
| exact, one orientation per call (today) | on | 12,298 | 12,326 | 4,500 | 402 | 43,416 | 98.731 |
| exact, 32 FZ points per call (BatchedScreen) | on | 12,298 | 12,326 | 4,500 | 402 | 43,416 | 98.731 |
| interpolated, one FZ point per call (CppOmegaInterp) | on | 12,325 | 12,326 | 4,885 | 1 | 43,970 | 99.991 |
| interpolated, 32 FZ points per call (both) | on | 12,325 | 12,326 | 4,885 | 1 | 43,969 | 99.989 |
<!-- /table:evalA_passsets_l0 -->

<!-- table:evalA_random_cpp -->
| voxel | find_py | find_cpp | cand_evaluated_py | conv_py | conv_cpp | var_py | var_cpp |
|---|---|---|---|---|---|---|---|
| 15901 | 193 | 580 | 1 | 1 | 1 | 1277 | 1894 |
| 17034 | 4376 | 162 | 30 | 0 | 1 | 1992 | 1388 |
| 18558 | 208 | 1528 | 1 | 1 | 1 | 1453 | 728 |
| 19889 | 161 | 162 | 1 | 1 | 1 | 1101 | 1245 |
| 22867 | 4205 | 162 | 30 | 0 | 1 | 1310 | 1564 |
| 5242 | 4984 | 162 | 30 | 0 | 1 | 1255 | 1718 |
<!-- /table:evalA_random_cpp -->

<!-- table:evalA_random_cpp_edge -->
| voxel | find_py | find_cpp | cand_evaluated_py | conv_py | conv_cpp | var_py | var_cpp |
|---|---|---|---|---|---|---|---|
| 15901 | 208 | 580 | 1 | 1 | 1 | 1948 | 1894 |
| 17034 | 161 | 162 | 1 | 1 | 1 | 1321 | 1388 |
| 18558 | 401 | 1528 | 2 | 1 | 1 | 1266 | 728 |
| 19889 | 161 | 162 | 1 | 1 | 1 | 1398 | 1245 |
| 22867 | 161 | 162 | 1 | 1 | 1 | 2025 | 1564 |
| 5242 | 97 | 162 | 1 | 1 | 1 | 771 | 1718 |
<!-- /table:evalA_random_cpp_edge -->

<!-- table:evalA_random_cpp_edge_interp -->
| voxel | find_py | find_cpp | cand_evaluated_py | conv_py | conv_cpp | var_py | var_cpp |
|---|---|---|---|---|---|---|---|
| 15901 | 193 | 580 | 1 | 1 | 1 | 2223 | 1894 |
| 17034 | 208 | 162 | 1 | 1 | 1 | 2245 | 1388 |
| 18558 | 208 | 1528 | 1 | 1 | 1 | 1255 | 728 |
| 19889 | 193 | 162 | 1 | 1 | 1 | 1002 | 1245 |
| 22867 | 161 | 162 | 1 | 1 | 1 | 2168 | 1564 |
| 5242 | 129 | 162 | 1 | 1 | 1 | 969 | 1718 |
<!-- /table:evalA_random_cpp_edge_interp -->

<!-- table:evalA_random_cpp_edge_interp_batched -->
| voxel | find_py | find_cpp | cand_evaluated_py | conv_py | conv_cpp | var_py | var_cpp |
|---|---|---|---|---|---|---|---|
| 15901 | 193 | 580 | 1 | 1 | 1 | 2223 | 1894 |
| 17034 | 208 | 162 | 1 | 1 | 1 | 2245 | 1388 |
| 18558 | 208 | 1528 | 1 | 1 | 1 | 1255 | 728 |
| 19889 | 193 | 162 | 1 | 1 | 1 | 1002 | 1245 |
| 22867 | 161 | 162 | 1 | 1 | 1 | 2168 | 1564 |
| 5242 | 129 | 162 | 1 | 1 | 1 | 969 | 1718 |
<!-- /table:evalA_random_cpp_edge_interp_batched -->

<!-- table:evalA_speed_screen -->
| arm | us_total | us_solve | us_overlap_c | us_rest | speedup |
|---|---|---|---|---|---|
| exact, one orientation per call (today) | 388.4 | 61.7 | 8.6 | 333.2 | 1.0 |
| exact, 32 FZ points per call (BatchedScreen) | 22.2 | 0.6 | 7.7 | 14.1 | 17.5 |
| interpolated, one FZ point per call (CppOmegaInterp) | 80.5 | 6.8 | 7.4 | 67.2 | 4.8 |
| interpolated, 32 FZ points per call (both) | 22.7 | 0.3 | 7.7 | 15.0 | 17.1 |
<!-- /table:evalA_speed_screen -->

<!-- table:evalA_speed_voxel -->
| switches | n | l0_screen_us | l0_screen_s | levels012_s | wall_s | evals | contended |
|---|---|---|---|---|---|---|---|
| off | 6 | 415.2 | 18.3 | 115.6 | 119.3 | 225,883 | 5 |
| BatchedScreen | 6 | 24.0 | 1.1 | 89.6 | 92.6 | 225,883 | 4 |
| CppOmegaInterp | 6 | 88.3 | 3.9 | 94.2 | 97.5 | 225,648 | 3 |
| both | 6 | 24.7 | 1.1 | 88.9 | 92.0 | 225,648 | 1 |
<!-- /table:evalA_speed_voxel -->

<!-- table:evalA_speed_voxel_quiet -->
| arm | voxel | label | wall_s | evals | l0_screen_us | us_per_eval | load |
|---|---|---|---|---|---|---|---|
| cpp_edge | 15901 | single-worker | 122.5 | 237,890 | 406.9 | 493 | 1.39 |
| cpp_edge | 5242 | single-worker | 100.8 | 199,386 | 398.7 | 486 | 1.29 |
| cpp_edge_interp_batched | 15901 | single-worker | 96.2 | 238,868 | 24.2 | 381 | 1.35 |
| cpp_edge_interp_batched | 5242 | single-worker | 75.6 | 199,195 | 24.1 | 360 | 1.19 |
<!-- /table:evalA_speed_voxel_quiet -->

<!-- table:evalA_stages_cpp -->
| stage | python | cpp | diff | share_of_gap |
|---|---|---|---|---|
| screen | 67,978 | 67,928 | 51 | -1 |
| recost | 26,148 | 26,046 | 102 | -1 |
| rerecost | 0 | 10,829 | -10,829 | 128 |
| quick | 130,266 | 129,948 | 318 | -4 |
| find | 2,354 | 459 | 1,895 | -22 |
| variance | 1,398 | 1,423 | -25 | 0 |
| final | 1 | 1 | 0 | -0 |
| total | 228,146 | 236,634 | -8,487 | 100 |
<!-- /table:evalA_stages_cpp -->

<!-- table:evalA_stages_cpp_edge -->
| stage | python | cpp | diff | share_of_gap |
|---|---|---|---|---|
| screen | 67,968 | 67,928 | 40 | -0 |
| recost | 26,077 | 26,046 | 31 | -0 |
| rerecost | 0 | 10,829 | -10,829 | 101 |
| quick | 130,184 | 129,948 | 236 | -2 |
| find | 198 | 459 | -261 | 2 |
| variance | 1,455 | 1,423 | 32 | -0 |
| final | 1 | 1 | 0 | -0 |
| total | 225,883 | 236,634 | -10,750 | 100 |
<!-- /table:evalA_stages_cpp_edge -->

<!-- table:evalA_stages_cpp_edge_interp -->
| stage | python | cpp | diff | share_of_gap |
|---|---|---|---|---|
| screen | 67,911 | 67,928 | -16 | 0 |
| recost | 26,016 | 26,046 | -30 | 0 |
| rerecost | 0 | 10,829 | -10,829 | 99 |
| quick | 129,894 | 129,948 | -54 | 0 |
| find | 182 | 459 | -277 | 3 |
| variance | 1,644 | 1,423 | 221 | -2 |
| final | 1 | 1 | 0 | -0 |
| total | 225,648 | 236,634 | -10,986 | 100 |
<!-- /table:evalA_stages_cpp_edge_interp -->

<!-- table:evalA_stages_cpp_edge_interp_batched -->
| stage | python | cpp | diff | share_of_gap |
|---|---|---|---|---|
| screen | 67,911 | 67,928 | -16 | 0 |
| recost | 26,016 | 26,046 | -30 | 0 |
| rerecost | 0 | 10,829 | -10,829 | 99 |
| quick | 129,894 | 129,948 | -54 | 0 |
| find | 182 | 459 | -277 | 3 |
| variance | 1,644 | 1,423 | 221 | -2 |
| final | 1 | 1 | 0 | -0 |
| total | 225,648 | 236,634 | -10,986 | 100 |
<!-- /table:evalA_stages_cpp_edge_interp_batched -->

<!-- table:evalA_units_cpp -->
| voxel | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | l0_py | l0_cpp | screen_per_fz_py | screen_per_fz_cpp | quick_per_cand_py | quick_per_cand_cpp | recost_per_pass_py | recost_per_pass_cpp | ratio_total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 7,708 | 7,726 | 28,920 | 28,873 | 11,498 | 11,562 | 8,110 | 8,151 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.946 |
| 17034 | 7,761 | 7,744 | 29,641 | 29,509 | 11,689 | 11,609 | 8,283 | 8,284 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.978 |
| 18558 | 7,344 | 7,336 | 21,922 | 21,860 | 10,004 | 9,973 | 7,048 | 7,026 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.954 |
| 19889 | 7,560 | 7,544 | 26,626 | 26,408 | 10,877 | 10,820 | 7,683 | 7,673 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.958 |
| 22867 | 7,739 | 7,724 | 29,663 | 29,472 | 11,613 | 11,549 | 8,233 | 8,234 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.973 |
| 5242 | 7,207 | 7,211 | 20,119 | 20,154 | 9,452 | 9,461 | 6,674 | 6,674 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.975 |
<!-- /table:evalA_units_cpp -->

<!-- table:evalA_units_cpp_edge -->
| voxel | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | l0_py | l0_cpp | screen_per_fz_py | screen_per_fz_cpp | quick_per_cand_py | quick_per_cand_cpp | recost_per_pass_py | recost_per_pass_cpp | ratio_total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 7,702 | 7,726 | 28,871 | 28,873 | 11,462 | 11,562 | 8,107 | 8,151 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.947 |
| 17034 | 7,751 | 7,744 | 29,632 | 29,509 | 11,646 | 11,609 | 8,264 | 8,284 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.956 |
| 18558 | 7,351 | 7,336 | 21,948 | 21,860 | 10,039 | 9,973 | 7,028 | 7,026 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.957 |
| 19889 | 7,556 | 7,544 | 26,527 | 26,408 | 10,857 | 10,820 | 7,688 | 7,673 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.958 |
| 22867 | 7,740 | 7,724 | 29,504 | 29,472 | 11,619 | 11,549 | 8,231 | 8,234 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.960 |
| 5242 | 7,212 | 7,211 | 19,981 | 20,154 | 9,469 | 9,461 | 6,668 | 6,674 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.950 |
<!-- /table:evalA_units_cpp_edge -->

<!-- table:evalA_units_cpp_edge_interp -->
| voxel | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | l0_py | l0_cpp | screen_per_fz_py | screen_per_fz_cpp | quick_per_cand_py | quick_per_cand_cpp | recost_per_pass_py | recost_per_pass_cpp | ratio_total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 7,712 | 7,726 | 28,959 | 28,873 | 11,507 | 11,562 | 8,149 | 8,151 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.951 |
| 17034 | 7,746 | 7,744 | 29,476 | 29,509 | 11,630 | 11,609 | 8,282 | 8,284 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.958 |
| 18558 | 7,353 | 7,336 | 21,849 | 21,860 | 10,040 | 9,973 | 7,026 | 7,026 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.955 |
| 19889 | 7,542 | 7,544 | 26,347 | 26,408 | 10,804 | 10,820 | 7,675 | 7,673 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.952 |
| 22867 | 7,717 | 7,724 | 29,438 | 29,472 | 11,530 | 11,549 | 8,231 | 8,234 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.955 |
| 5242 | 7,204 | 7,211 | 20,028 | 20,154 | 9,436 | 9,461 | 6,674 | 6,674 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.949 |
<!-- /table:evalA_units_cpp_edge_interp -->

<!-- table:evalA_units_cpp_edge_interp_batched -->
| voxel | fz_py | fz_cpp | passed_py | passed_cpp | cands_py | cands_cpp | l0_py | l0_cpp | screen_per_fz_py | screen_per_fz_cpp | quick_per_cand_py | quick_per_cand_cpp | recost_per_pass_py | recost_per_pass_cpp | ratio_total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15901 | 7,712 | 7,726 | 28,959 | 28,873 | 11,507 | 11,562 | 8,149 | 8,151 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.951 |
| 17034 | 7,746 | 7,744 | 29,476 | 29,509 | 11,630 | 11,609 | 8,282 | 8,284 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.958 |
| 18558 | 7,353 | 7,336 | 21,849 | 21,860 | 10,040 | 9,973 | 7,026 | 7,026 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.955 |
| 19889 | 7,542 | 7,544 | 26,347 | 26,408 | 10,804 | 10,820 | 7,675 | 7,673 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.952 |
| 22867 | 7,717 | 7,724 | 29,438 | 29,472 | 11,530 | 11,549 | 8,231 | 8,234 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.955 |
| 5242 | 7,204 | 7,211 | 20,028 | 20,154 | 9,436 | 9,461 | 6,674 | 6,674 | 9.00 | 9.00 | 12.00 | 12.00 | 1.00 | 1.00 | 0.949 |
<!-- /table:evalA_units_cpp_edge_interp_batched -->
