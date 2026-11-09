<!-- table:cross_old_to_cpp -->
| level | c | n_stop | broken_1deg | broken_2deg | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|
| 0 | 0.75 | 14 | 0 | 0 | 0.27 | 0.16 |
| 1 | 0.75 | 64 | 0 | 0 | 0.62 | 0.17 |
| 2 | 0.70 | 108 | 0 | 0 | 0.83 | 0.06 |
<!-- /table:cross_old_to_cpp -->

<!-- table:final_sep_cpp -->
| c | right_le | wrong_le | right_gt | wrong_gt |
|---|---|---|---|---|
| 0.20 | 130 | 0 | 14 | 56 |
| 0.30 | 141 | 0 | 3 | 56 |
| 0.40 | 144 | 0 | 0 | 56 |
| 0.50 | 144 | 0 | 0 | 56 |
| 0.60 | 144 | 0 | 0 | 56 |
| 0.65 | 144 | 0 | 0 | 56 |
| 0.70 | 144 | 8 | 0 | 48 |
| 0.75 | 144 | 28 | 0 | 28 |
| 0.80 | 144 | 40 | 0 | 16 |
| 0.85 | 144 | 48 | 0 | 8 |
| 0.90 | 144 | 51 | 0 | 5 |
| 0.95 | 144 | 56 | 0 | 0 |
<!-- /table:final_sep_cpp -->

<!-- table:final_sep_old -->
| c | right_le | wrong_le | right_gt | wrong_gt |
|---|---|---|---|---|
| 0.20 | 149 | 0 | 10 | 41 |
| 0.30 | 153 | 0 | 6 | 41 |
| 0.40 | 159 | 0 | 0 | 41 |
| 0.50 | 159 | 0 | 0 | 41 |
| 0.60 | 159 | 0 | 0 | 41 |
| 0.65 | 159 | 0 | 0 | 41 |
| 0.70 | 159 | 4 | 0 | 37 |
| 0.75 | 159 | 24 | 0 | 17 |
| 0.80 | 159 | 34 | 0 | 7 |
| 0.85 | 159 | 36 | 0 | 5 |
| 0.90 | 159 | 39 | 0 | 2 |
| 0.95 | 159 | 41 | 0 | 0 |
<!-- /table:final_sep_old -->

<!-- table:levels_cpp -->
| level | best_le1_right | best_le2_right | best_le2_wrong | kept_le2_right | kept_le2_wrong | far_cost_min | evals_after_pct | evals_after_right_mean |
|---|---|---|---|---|---|---|---|---|
| 0 | 69 | 72 | 1 | 137 | 14 | 0.772 | 3.4 | 1686 |
| 1 | 119 | 119 | 1 | 144 | 4 | 0.757 | 1.9 | 938 |
| 2 | 139 | 139 | 0 | 144 | 0 | 0.735 | 1.5 | 733 |
| 3 | 143 | 143 | 0 | 143 | 0 | 0.735 | 1.3 | 683 |
<!-- /table:levels_cpp -->

<!-- table:levels_old -->
| level | best_le1_right | best_le2_right | best_le2_wrong | kept_le2_right | kept_le2_wrong | far_cost_min | evals_after_pct | evals_after_right_mean |
|---|---|---|---|---|---|---|---|---|
| 0 | 71 | 74 | 1 | 156 | 14 | 0.777 | 4.8 | 2411 |
| 1 | 141 | 141 | 1 | 159 | 7 | 0.773 | 2.6 | 1292 |
| 2 | 148 | 148 | 0 | 159 | 1 | 0.690 | 1.9 | 967 |
| 3 | 155 | 155 | 0 | 155 | 0 | 0.690 | 1.7 | 886 |
<!-- /table:levels_old -->

<!-- table:pilot_hit -->
| run | group | n | n_wrong | hit_right_min | hit_right_median | hit_wrong_median | hit_wrong_max | evals_in_wrong_pct |
|---|---|---|---|---|---|---|---|---|
| keyoff_rerun | full_search | 69 | 4 | 0.936 | 1.000 | 0.777 | 0.904 | 5.8 |
| keyoff_rerun | local_fit | 1931 | 33 | 0.888 | 0.994 | 0.904 | 0.918 | 1.2 |
| outward_only | full_search | 70 | 6 | 0.951 | 1.000 | 0.875 | 0.892 | 2.1 |
| outward_only | local_fit | 1930 | 30 | 0.900 | 0.992 | 0.902 | 0.918 | 0.9 |
| integrated | full_search | 70 | 4 | 0.964 | 1.000 | 0.879 | 0.899 | 0.4 |
| integrated | local_fit | 1930 | 38 | 0.900 | 0.994 | 0.902 | 0.918 | 1.4 |
| smooth_only | full_search | 65 | 0 | 0.921 | 1.000 | n/a | n/a | 0.0 |
| smooth_only | local_fit | 1935 | 34 | 0.859 | 0.992 | 0.902 | 0.918 | 1.2 |
<!-- /table:pilot_hit -->

<!-- table:pilot_thr_full_search -->
| h | right_ge | wrong_ge | right_lt | wrong_lt |
|---|---|---|---|---|
| 0.90 | 65 | 1 | 0 | 3 |
| 0.93 | 65 | 0 | 0 | 4 |
| 0.95 | 64 | 0 | 1 | 4 |
| 0.96 | 62 | 0 | 3 | 4 |
| 0.97 | 59 | 0 | 6 | 4 |
| 0.98 | 57 | 0 | 8 | 4 |
| 0.99 | 44 | 0 | 21 | 4 |
<!-- /table:pilot_thr_full_search -->

<!-- table:pilot_thr_local_fit -->
| h | right_ge | wrong_ge | right_lt | wrong_lt |
|---|---|---|---|---|
| 0.90 | 1896 | 28 | 2 | 5 |
| 0.93 | 1883 | 0 | 15 | 33 |
| 0.95 | 1852 | 0 | 46 | 33 |
| 0.96 | 1789 | 0 | 109 | 33 |
| 0.97 | 1590 | 0 | 308 | 33 |
| 0.98 | 1393 | 0 | 505 | 33 |
| 0.99 | 1076 | 0 | 822 | 33 |
<!-- /table:pilot_thr_local_fit -->

<!-- table:safe_cpp -->
| max_broken | level | c | n_stop | broken_1deg | broken_2deg | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0.75 | 14 | 0 | 0 | 0.27 | 0.16 |
| 0 | 1 | 0.75 | 64 | 0 | 0 | 0.62 | 0.17 |
| 0 | 2 | 0.75 | 117 | 0 | 0 | 0.90 | 0.06 |
| 1 | 0 | 0.80 | 19 | 1 | 1 | 0.36 | 0.22 |
| 1 | 1 | 0.75 | 64 | 0 | 0 | 0.62 | 0.17 |
| 1 | 2 | 0.80 | 135 | 1 | 1 | 1.04 | 0.07 |
<!-- /table:safe_cpp -->

<!-- table:safe_old -->
| max_broken | level | c | n_stop | broken_1deg | broken_2deg | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0.75 | 17 | 0 | 0 | 0.39 | 0.25 |
| 0 | 1 | 0.75 | 88 | 0 | 0 | 1.16 | 0.37 |
| 0 | 2 | 0.70 | 127 | 0 | 0 | 1.26 | 0.11 |
| 1 | 0 | 0.80 | 28 | 1 | 1 | 0.70 | 0.44 |
| 1 | 1 | 0.75 | 88 | 0 | 0 | 1.16 | 0.37 |
| 1 | 2 | 0.70 | 127 | 0 | 0 | 1.26 | 0.11 |
<!-- /table:safe_old -->

<!-- table:sep_levels_cpp -->
| level | right_max | wrong_min | n_right | n_wrong | separable |
|---|---|---|---|---|---|
| L0 | 0.963 | 0.792 | 144 | 56 | False |
| L1 | 0.940 | 0.757 | 144 | 56 | False |
| L2 | 0.917 | 0.735 | 144 | 56 | False |
| L3 | 0.867 | 0.735 | 144 | 56 | False |
<!-- /table:sep_levels_cpp -->

<!-- table:sep_levels_old -->
| level | right_max | wrong_min | n_right | n_wrong | separable |
|---|---|---|---|---|---|
| L0 | 0.957 | 0.850 | 159 | 41 | False |
| L1 | 0.938 | 0.805 | 159 | 41 | False |
| L2 | 0.891 | 0.690 | 159 | 41 | False |
| L3 | 0.879 | 0.690 | 159 | 41 | False |
<!-- /table:sep_levels_old -->

<!-- table:split_after -->
| after_level | e0_evals_pct | full_evals_pct | full_time_pct |
|---|---|---|---|
| 0 | 4.8 | 39.2 | 38.2 |
| 1 | 2.6 | 12.9 | 12.7 |
| 2 | 1.9 | 4.1 | 4.3 |
| 3 | 1.7 | 1.3 | 1.6 |
<!-- /table:split_after -->

<!-- table:split_stages -->
| stage | e0_evals_pct | full_evals_pct | full_time_pct |
|---|---|---|---|
| level 0 screen and re-cost | 90.7 | 23.7 | 15.8 |
| level 0 quick MC | 4.5 | 37.1 | 46.0 |
| level 1 | 2.3 | 26.3 | 25.4 |
| level 2 | 0.7 | 8.8 | 8.4 |
| level 3 | 0.2 | 2.8 | 2.7 |
| finisher | 1.7 | 1.3 | 1.6 |
| not in a stage (time only) | n/a | n/a | 0.0 |
<!-- /table:split_stages -->

<!-- table:stages_cpp -->
| level | diameter_deg | screen_pct | quick_pct | n_disc_median | right_with_basin | right_no_basin | wrong_with_basin | basin_pos_median | wrong_basin_pos_median | basin_frac_median | basin_frac_p90 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 5.00 | 93.4 | 3.2 | 125 | 144 | 0 | 56 | 5 | 10 | 0.04 | 0.18 |
| 1 | 3.33 | 0.7 | 0.8 | 32 | 144 | 0 | 38 | 0 | 7 | 0.00 | 0.21 |
| 2 | 2.22 | 0.2 | 0.2 | 8 | 144 | 0 | 5 | 0 | 3 | 0.00 | 0.12 |
| 3 | 1.48 | 0.1 | 0.0 | 2 | 144 | 0 | 0 | 0 | n/a | 0.00 | 0.00 |
<!-- /table:stages_cpp -->

<!-- table:stages_old -->
| level | diameter_deg | screen_pct | quick_pct | n_disc_median | right_with_basin | right_no_basin | wrong_with_basin | basin_pos_median | wrong_basin_pos_median | basin_frac_median | basin_frac_p90 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 5.00 | 90.7 | 4.5 | 182 | 159 | 0 | 41 | 3 | 6 | 0.02 | 0.10 |
| 1 | 3.33 | 1.1 | 1.2 | 48 | 159 | 0 | 28 | 0 | 8 | 0.00 | 0.05 |
| 2 | 2.22 | 0.4 | 0.3 | 12 | 159 | 0 | 7 | 0 | 3 | 0.00 | 0.00 |
| 3 | 1.48 | 0.1 | 0.1 | 3 | 159 | 0 | 1 | 0 | 2 | 0.00 | 0.00 |
<!-- /table:stages_old -->

<!-- table:stop_grid_cpp -->
| level | c | n_stop | broken_1deg | broken_2deg | stop_wrong_base | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|---|
| 0 | 0.30 | 2 | 0 | 0 | 0 | 0.03 | 0.02 |
| 0 | 0.40 | 2 | 0 | 0 | 0 | 0.03 | 0.02 |
| 0 | 0.50 | 2 | 0 | 0 | 0 | 0.03 | 0.02 |
| 0 | 0.60 | 5 | 0 | 0 | 0 | 0.09 | 0.06 |
| 0 | 0.70 | 9 | 0 | 0 | 0 | 0.17 | 0.10 |
| 0 | 0.80 | 19 | 1 | 1 | 1 | 0.36 | 0.22 |
| 0 | 0.90 | 70 | 16 | 16 | 6 | 1.25 | 0.75 |
| 1 | 0.30 | 3 | 0 | 0 | 0 | 0.03 | 0.01 |
| 1 | 0.40 | 5 | 0 | 0 | 0 | 0.04 | 0.01 |
| 1 | 0.50 | 12 | 0 | 0 | 0 | 0.11 | 0.03 |
| 1 | 0.60 | 30 | 0 | 0 | 0 | 0.29 | 0.08 |
| 1 | 0.70 | 51 | 0 | 0 | 0 | 0.50 | 0.14 |
| 1 | 0.80 | 84 | 3 | 3 | 3 | 0.82 | 0.22 |
| 1 | 0.90 | 147 | 19 | 19 | 18 | 1.43 | 0.39 |
| 2 | 0.30 | 7 | 0 | 0 | 0 | 0.05 | 0.00 |
| 2 | 0.40 | 22 | 0 | 0 | 0 | 0.16 | 0.01 |
| 2 | 0.50 | 46 | 0 | 0 | 0 | 0.35 | 0.02 |
| 2 | 0.60 | 77 | 0 | 0 | 0 | 0.58 | 0.04 |
| 2 | 0.70 | 108 | 0 | 0 | 0 | 0.83 | 0.06 |
| 2 | 0.80 | 135 | 1 | 1 | 7 | 1.04 | 0.07 |
| 2 | 0.90 | 173 | 4 | 4 | 32 | 1.29 | 0.09 |
<!-- /table:stop_grid_cpp -->

<!-- table:stop_grid_old -->
| level | c | n_stop | broken_1deg | broken_2deg | stop_wrong_base | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|---|
| 0 | 0.30 | 2 | 0 | 0 | 0 | 0.05 | 0.03 |
| 0 | 0.40 | 2 | 0 | 0 | 0 | 0.05 | 0.03 |
| 0 | 0.50 | 2 | 0 | 0 | 0 | 0.05 | 0.03 |
| 0 | 0.60 | 5 | 0 | 0 | 0 | 0.12 | 0.08 |
| 0 | 0.70 | 9 | 0 | 0 | 0 | 0.21 | 0.13 |
| 0 | 0.80 | 28 | 1 | 1 | 0 | 0.70 | 0.44 |
| 0 | 0.90 | 83 | 27 | 27 | 3 | 2.06 | 1.30 |
| 1 | 0.30 | 6 | 0 | 0 | 0 | 0.08 | 0.02 |
| 1 | 0.40 | 7 | 0 | 0 | 0 | 0.09 | 0.03 |
| 1 | 0.50 | 19 | 0 | 0 | 0 | 0.25 | 0.08 |
| 1 | 0.60 | 38 | 0 | 0 | 0 | 0.51 | 0.16 |
| 1 | 0.70 | 61 | 0 | 0 | 0 | 0.80 | 0.25 |
| 1 | 0.80 | 107 | 2 | 2 | 0 | 1.42 | 0.45 |
| 1 | 0.90 | 170 | 16 | 16 | 19 | 2.22 | 0.70 |
| 2 | 0.30 | 16 | 0 | 0 | 0 | 0.15 | 0.01 |
| 2 | 0.40 | 30 | 0 | 0 | 0 | 0.29 | 0.02 |
| 2 | 0.50 | 55 | 0 | 0 | 0 | 0.54 | 0.05 |
| 2 | 0.60 | 85 | 0 | 0 | 0 | 0.84 | 0.07 |
| 2 | 0.70 | 127 | 0 | 0 | 1 | 1.26 | 0.11 |
| 2 | 0.80 | 155 | 7 | 7 | 6 | 1.54 | 0.13 |
| 2 | 0.90 | 192 | 11 | 11 | 33 | 1.85 | 0.16 |
<!-- /table:stop_grid_old -->

<!-- table:zero_counts -->
| count | k | n | lo_pct | hi_pct |
|---|---|---|---|---|
| E0 old rule: voxels misclassified by the final cost (best c) | 0 | 200 | 0.0 | 1.9 |
| E0 C++ rule: voxels misclassified by the final cost (best c) | 0 | 200 | 0.0 | 1.9 |
| E0 old rule: broken at 1 degree, safe rule after level 0 | 0 | 159 | 0.0 | 2.4 |
| E0 old rule: broken at 1 degree, safe rule after level 1 | 0 | 159 | 0.0 | 2.4 |
| E0 old rule: broken at 1 degree, safe rule after level 2 | 0 | 159 | 0.0 | 2.4 |
| E0 C++ rule: broken at 1 degree, old-rule c after level 0 | 0 | 144 | 0.0 | 2.6 |
| E0 C++ rule: broken at 1 degree, old-rule c after level 1 | 0 | 144 | 0.0 | 2.6 |
| E0 C++ rule: broken at 1 degree, old-rule c after level 2 | 0 | 144 | 0.0 | 2.6 |
| pilot key-off: wrong voxels at or above the cut | 0 | 37 | 0.0 | 9.4 |
| pilot key-off: answers at or above the cut that are wrong | 0 | 1948 | 0.0 | 0.2 |
<!-- /table:zero_counts -->
