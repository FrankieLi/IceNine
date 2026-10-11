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
